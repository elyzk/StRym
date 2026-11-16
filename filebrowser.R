# filebrowser.R -- StRym's own file and folder chooser (the "Browse…" window), the same on Windows, macOS and Linux.
#
# Why not the usual shinyFiles window: it asks the operating system about every drive each time it opens (PowerShell on Windows takes seconds) and it
# looks inside every sub-folder it shows, which is very slow on external drives and on Windows folders seen from WSL. This chooser
#   * opens at once: the list of places on the left (Desktop, Documents, Downloads, every drive, the whole file system, folders you used before) is built from
#     what is already known, and the answer from Windows about unplugged or unmounted drives arrives a moment later without making you wait;
#   * lists a folder with ONE quick call and shows sizes only for the first few hundred files, so a folder with thousands of files appears immediately;
#   * has a filter box (it filters in your browser, no waiting), a "Show all files" switch (the list normally shows only the file types the box expects),
#     an Up button and a path box (type or paste a path, press Enter);
#   * mounts a Windows drive for you when you click it and it is not connected yet (see drives.R).
#
# Use:  br <- strym_browser_attach(input, output, session, key = "bf")        once per Shiny module / session
#       br$open(mode = "file" | "dir" | "any", title, exts, multiple, start, cb = function(paths) ...)
# `cb` receives the chosen file path(s) or folder. pathUI()/pathServer() (gui_helpers.R) and the TUX tab (app.R) use it.

.bf_recent <- new.env()                        # folders used before in this app session (newest first), shared by every Browse window
.bf_cache <- new.env()                         # small caches (the Windows profile folder)
BF_MAX_ROWS <- 1500                            # rows drawn per folder; the filter box searches only these, the note says how many there were
BF_STAT_MAX <- 300                             # sizes and dates are read for at most this many files (each costs a call on a slow drive)
BF_HIDE <- c("$RECYCLE.BIN", "System Volume Information", "desktop.ini", "Thumbs.db", "pagefile.sys", "hiberfil.sys", "swapfile.sys", "DumpStack.log.tmp", "lost+found")

# The folder with the Windows user's Desktop / Documents / Downloads when StRym runs inside WSL (best guess: the profile used most recently).
.win_profile <- function() {
  if (!is.null(.bf_cache$prof) && as.numeric(Sys.time()) - .bf_cache$prof$t < 120) return(.bf_cache$prof$v)
  v <- NULL; base <- "/mnt/c/Users"
  if (dir.exists(base)) {
    u <- list.dirs(base, recursive = FALSE, full.names = TRUE); u <- u[!basename(u) %in% c("Public", "Default", "Default User", "All Users")]
    u <- u[dir.exists(file.path(u, "Desktop")) | dir.exists(file.path(u, "Documents")) | dir.exists(file.path(u, "Downloads"))]
    if (length(u)) v <- u[which.max(suppressWarnings(file.mtime(u)))]
  }
  .bf_cache$prof <- list(t = as.numeric(Sys.time()), v = v); v
}

# The places shown on the left. Returns data.frame(group, label, path, kind, letter, note); kind "go" opens the folder, "mount" connects a Windows drive first.
strym_places <- function(drives = drives_nowait()) {
  pl <- strym_platform(); rows <- list()
  add <- function(group, label, path, kind = "go", letter = "", note = "") rows[[length(rows) + 1]] <<- data.frame(group = group, label = label, path = path, kind = kind,
                                                                                                         letter = letter, note = note, stringsAsFactors = FALSE)
  home <- path.expand("~"); prof <- if (pl == "wsl") .win_profile() else home
  for (nm in c("Desktop", "Documents", "Downloads")) { p <- if (!is.null(prof)) file.path(prof, nm) else ""; if (nzchar(p) && dir.exists(p)) add("Quick places", nm, p) }
  if (pl == "wsl") { if (!is.null(prof)) add("Quick places", "Windows user folder", prof); add("Quick places", "Linux home folder (StRym)", home)
  } else add("Quick places", "Home", home)
  data_root <- Sys.getenv("STRYM_DATA_ROOT", "")
  if (nzchar(data_root) && dir.exists(data_root)) add("Quick places", "Data", path.expand(data_root))
  if (is.data.frame(drives) && nrow(drives)) {
    drives <- drives[drives$state != "removed", , drop = FALSE]
    for (i in seq_len(nrow(drives))) {
      d <- drives[i, ]; free <- if (is.finite(d$free)) paste(fmt_bytes(d$free), "free") else ""
      lab <- if (pl == "wsl") sprintf("%s:  %s", d$letter, if (nzchar(d$label)) d$label else if (d$letter == "C") "Windows" else "drive") else if (nzchar(d$label)) d$label else basename(d$mount)
      if (d$state == "ready" || pl != "wsl") add("Drives", lab, d$mount, "go", d$letter, if (d$state == "ready") free else "not responding")
      else add("Drives", lab, d$mount, "mount", d$letter, if (d$state == "stale") "click to reconnect" else "click to connect")
    }
  }
  add("Drives", switch(pl, mac = "Macintosh HD (whole computer)", windows = "C:  (whole drive)", "Whole file system  /"), if (pl == "windows") "C:/" else "/")
  r <- do.call(rbind, rows)
  r[!duplicated(r$path) | r$kind == "mount", , drop = FALSE]
}
.bf_remember <- function(dir) {                                       # newest first; a folder used again moves to the top
  if (!length(dir) || !nzchar(dir[1]) || !dir.exists(dir[1])) return(invisible())
  k <- ls(.bf_recent); if (length(k)) { v <- unlist(mget(k, envir = .bf_recent), use.names = FALSE); rm(list = k[v == dir[1]], envir = .bf_recent) }
  assign(sprintf("%.3f", as.numeric(Sys.time())), dir[1], envir = .bf_recent)
}
.bf_recent_dirs <- function(n = 5) { k <- ls(.bf_recent); if (!length(k)) return(character()); v <- unlist(mget(k, envir = .bf_recent), use.names = FALSE); utils::head(v[order(-as.numeric(k))], n) }

# One quick listing of a folder: list(ok, error, entries = data.frame(name, dir, size, mtime), total, shown).
# `exts` limits the FILES shown (NULL = all); folders are always shown. Uses a single `ls` call (a time limit stops a drive that does not answer).
strym_list_dir <- function(path, exts = NULL, show_all = FALSE, hidden = FALSE, timeout = 20, max_rows = BF_MAX_ROWS) {
  bad <- function(msg) list(ok = FALSE, error = msg, entries = NULL, total = 0L, shown = 0L)
  if (!nzchar(path) || !dir.exists(path)) return(bad(if (nzchar(path)) "This folder does not exist (is the drive connected?)." else "Choose a place on the left."))
  # ONE call lists the folder. GNU find (Linux, WSL) reads the entry types straight from the directory (a folder with 5,000 files in 0.2 s; `ls -p` needs a call
  # per file and took 4 s on a Windows folder seen from WSL). macOS has no GNU find, so it uses `ls -p` (a Mac reads its disks quickly).
  gnu <- !identical(strym_platform(), "mac")
  r <- tryCatch(if (gnu) processx::run("find", c(path, "-mindepth", "1", "-maxdepth", "1", "-printf", "%y\\t%f\\n"), error_on_status = FALSE, timeout = timeout)
                else processx::run("ls", c("-1Ap", "--", path), error_on_status = FALSE, timeout = timeout), error = function(e) NULL)
  if (is.null(r) || isTRUE(r$timeout)) return(bad(sprintf("This folder did not answer within %d seconds. If it is on an external drive, check that the drive is plugged in and not asleep, then press the refresh button.", timeout)))
  out <- enc2utf8(strsplit(r$stdout, "\n", fixed = TRUE)[[1]]); out <- out[nzchar(out)]
  if (!length(out) && !identical(r$status, 0L)) return(bad(if (grepl("ermission", r$stderr)) "You do not have permission to open this folder." else "This folder cannot be read."))
  if (gnu) { ty <- substr(out, 1, 1); nm <- substring(out, 3); isdir <- ty == "d"
             lk <- which(ty == "l"); if (length(lk) && length(lk) <= 300) isdir[lk] <- dir.exists(file.path(path, nm[lk]))        # a link to a folder is a folder
  } else { isdir <- grepl("/$", out); nm <- sub("/$", "", out) }
  keep <- !nm %in% BF_HIDE & (hidden | !startsWith(nm, "."))
  if (!is.null(exts) && !show_all) {
    ex <- tolower(tools::file_ext(nm)); ex2 <- tolower(tools::file_ext(sub("\\.gz$", "", nm, ignore.case = TRUE)))
    keep <- keep & (isdir | ex %in% tolower(exts) | ex2 %in% tolower(exts))
  }
  nm <- nm[keep]; isdir <- isdir[keep]
  o <- order(!isdir, tolower(nm), method = "radix"); nm <- nm[o]; isdir <- isdir[o]
  total <- length(nm); if (total > max_rows) { nm <- nm[seq_len(max_rows)]; isdir <- isdir[seq_len(max_rows)] }
  size <- rep(NA_real_, length(nm)); mt <- rep(NA_real_, length(nm))
  fi <- which(!isdir); if (length(fi) && length(fi) <= BF_STAT_MAX) { inf <- suppressWarnings(file.info(file.path(path, nm[fi]), extra_cols = FALSE)); size[fi] <- inf$size; mt[fi] <- as.numeric(inf$mtime) }
  list(ok = TRUE, error = "", entries = data.frame(name = nm, dir = isdir, size = size, mtime = mt, stringsAsFactors = FALSE), total = total, shown = length(nm))
}

# options(strym.debug = TRUE) (or the environment variable STRYM_DEBUG=1) prints a time stamp at the main steps of opening a window: where does the time go?
.bf_t <- function(what) if (isTRUE(getOption("strym.debug")) || nzchar(Sys.getenv("STRYM_DEBUG"))) message(sprintf("[browse %s] %s", format(Sys.time(), "%H:%M:%OS3"), what))

# ---- the window ----------------------------------------------------------------------------------------------------------------------------------
strym_browser_attach <- function(input, output, session, key = "bf") {
  ns <- session$ns; id <- function(x) paste0(key, "_", x)
  st <- reactiveValues(open = FALSE, path = "", mode = "dir", exts = NULL, multiple = FALSE, show_all = FALSE, hidden = FALSE, title = "", tick = 0, msg = "")
  cb <- NULL; sel <- reactiveVal(character())
  go_to <- function(p) { p <- sub("(?<=.)/+$", "", p, perl = TRUE); if (dir.exists(p)) { st$path <- normalizePath(p, winslash = "/", mustWork = FALSE); st$msg <- ""; sel(character())
                                                                              updateTextInput(session, id("path"), value = if (strym_is_wsl()) to_windows_path(st$path) else st$path) } else st$msg <- "That folder does not exist." }
  listing <- reactive({ req(st$open, nzchar(st$path)); st$tick; strym_list_dir(st$path, if (st$mode == "dir") NULL else st$exts, st$show_all, st$hidden) })

  output[[id("places")]] <- renderUI({
    req(st$open); invalidateLater(3000)                                  # the answer about drives that are not connected yet arrives in the background: pick it up
    .bf_t("places: start"); pl <- strym_places(); rec <- .bf_recent_dirs(); .bf_t("places: computed")
    groups <- c("Quick places", "Drives"); blocks <- lapply(groups, function(g) { x <- pl[pl$group == g, , drop = FALSE]; if (!nrow(x)) return(NULL)
      tagList(div(class = "bf-gh", g), lapply(seq_len(nrow(x)), function(i) div(class = paste("bf-place", if (identical(x$kind[i], "mount")) "bf-mount"), `data-input` = ns(id("pl")),
        `data-kind` = x$kind[i], `data-path` = x$path[i], `data-letter` = x$letter[i], title = x$path[i],
        span(class = "bf-pname", x$label[i]), if (nzchar(x$note[i])) span(class = "bf-pnote", x$note[i])))) })
    recent <- if (length(rec)) tagList(div(class = "bf-gh", "Recent folders"), lapply(rec, function(p) div(class = "bf-place", `data-input` = ns(id("pl")), `data-kind` = "go", `data-path` = p,
                                                                                                            `data-letter` = "", title = p, span(class = "bf-pname", basename(p)), span(class = "bf-pnote", dirname(p)))))
    div(class = "bf-places", blocks, recent, div(class = "bf-place bf-refresh", `data-input` = ns(id("pl")), `data-kind` = "refresh", `data-path` = "", `data-letter` = "", span(class = "bf-pname", "↻ Look for drives again")))
  })
  output[[id("list")]] <- renderUI({
    .bf_t("list: start"); L <- listing(); .bf_t("list: read"); if (!L$ok) return(div(class = "bf-empty bad", L$error))
    e <- L$entries; if (!nrow(e)) return(div(class = "bf-empty", if (st$mode == "dir") "This folder has no sub-folders." else "No matching files here. Tick “Show all files” to see everything."))
    pickable <- if (st$mode == "dir") rep(FALSE, nrow(e)) else !e$dir
    rows <- lapply(seq_len(nrow(e)), function(i) div(class = paste("bf-row", if (e$dir[i]) "bf-dir" else if (pickable[i]) "bf-file" else "bf-dim"), `data-name` = e$name[i],
      span(class = "bf-ico", if (e$dir[i]) "📁" else "📄"), span(class = "bf-name", e$name[i]),
      span(class = "bf-size", if (is.finite(e$size[i])) fmt_bytes(e$size[i]) else ""), span(class = "bf-time", if (is.finite(e$mtime[i])) format(as.POSIXct(e$mtime[i], origin = "1970-01-01"), "%Y-%m-%d %H:%M") else "")))
    div(class = "bf-list", `data-cd` = ns(id("cd")), `data-sel` = ns(id("sel")), `data-multi` = tolower(as.character(st$multiple)), `data-base` = st$path, rows)
  })
  output[[id("status")]] <- renderUI({
    L <- listing(); n <- length(sel())
    tagList(if (nzchar(st$msg)) div(class = "bad", st$msg),
            if (L$ok && L$total > L$shown) div(class = "hint", sprintf("%s items here; the first %s are listed. Use the path box to go deeper.", format(L$total, big.mark = ","), format(L$shown, big.mark = ","))),
            if (st$mode != "dir") div(class = "hint", if (n) sprintf("%d file%s selected.", n, if (n == 1) "" else "s") else if (st$multiple) "Click files to select them (Ctrl or Shift for several)." else "Click a file to select it."))
  })

  observeEvent(input[[id("cd")]], go_to(input[[id("cd")]]))
  observeEvent(input[[id("up")]], go_to(dirname(st$path)))
  observeEvent(input[[id("go")]], { p <- clean_path(input[[id("path")]]); if (dir.exists(p)) go_to(p) else if (st$mode != "dir" && file.exists(p)) { sel(basename(p)); go_to(dirname(p)); st$msg <- "" } else st$msg <- "That folder does not exist." })
  observeEvent(input[[id("refresh")]], { st$tick <- st$tick + 1; refresh_drives() })
  observeEvent(input[[id("all")]], st$show_all <- isTRUE(input[[id("all")]]), ignoreInit = TRUE)
  observeEvent(input[[id("hid")]], st$hidden <- isTRUE(input[[id("hid")]]), ignoreInit = TRUE)
  observeEvent(input[[id("sel")]], sel(as.character(unlist(input[[id("sel")]]))), ignoreNULL = FALSE)
  observeEvent(input[[id("pl")]], {
    x <- input[[id("pl")]]; k <- x$kind
    if (identical(k, "refresh")) { refresh_drives(); st$tick <- st$tick + 1; showNotification("Looking for drives…", duration = 3); return() }
    if (identical(k, "mount")) {
      nid <- showNotification(sprintf("Connecting %s: …", x$letter), duration = NULL, closeButton = FALSE, type = "message")
      r <- if (strym_is_wsl()) mount_drive(x$letter) else list(ok = dir.exists(x$path), message = "This drive is not available.")
      removeNotification(nid); if (isTRUE(r$ok)) { go_to(x$path) } else st$msg <- r$message; return()
    }
    go_to(x$path)
  })
  pick <- function(paths) { paths <- paths[nzchar(paths)]; if (!length(paths)) return(); .bf_remember(if (st$mode == "dir" || dir.exists(paths[1])) paths[1] else dirname(paths[1]))
                            removeModal(session); st$open <- FALSE; if (is.function(cb)) cb(paths) }
  observeEvent(input[[id("pick")]], { if (length(sel())) pick(file.path(st$path, sel())) })
  observeEvent(input[[id("pickdir")]], { p <- clean_path(input[[id("path")]]); if (!nzchar(p)) p <- st$path; pick(p) })
  observeEvent(input[[id("cancel")]], { removeModal(session); st$open <- FALSE })
  observe({ session$sendCustomMessage("tux-toggle", list(id = ns(id("pick")), enable = length(sel()) > 0)) })

  open <- function(mode = c("file", "dir", "any"), title = "Choose", exts = NULL, multiple = FALSE, start = "", cb = NULL) {
    .bf_t("open: called"); mode <- match.arg(mode); cb <<- cb
    p <- clean_path(start); if (nzchar(p) && !dir.exists(p)) p <- if (nzchar(dirname(p)) && dir.exists(dirname(p))) dirname(p) else ""
    if (!nzchar(p)) { r <- .bf_recent_dirs(1); p <- if (length(r)) r else { pl <- strym_places(); dl <- pl$path[pl$label == "Downloads"]; if (length(dl)) dl[1] else path.expand("~") } }
    st$mode <- mode; st$exts <- exts; st$multiple <- multiple; st$show_all <- FALSE; st$title <- title; st$msg <- ""; st$open <- TRUE; sel(character())
    st$path <- normalizePath(p, winslash = "/", mustWork = FALSE)
    showModal(modalDialog(title = title, size = "l", easyClose = FALSE, class = "bf-modal",
      div(class = "bf",
        div(class = "bf-side", uiOutput(ns(id("places")))),
        div(class = "bf-main",
          div(class = "bf-bar", actionButton(ns(id("up")), "↑ Up", class = "btn-default btn-sm"),
              textInput(ns(id("path")), NULL, value = if (strym_is_wsl()) to_windows_path(st$path) else st$path, width = "100%"),
              actionButton(ns(id("go")), "Go", class = "btn-default btn-sm"), actionButton(ns(id("refresh")), "↻", class = "btn-default btn-sm", title = "Read this folder again")),
          div(class = "bf-tools", tags$input(type = "search", class = "bf-filter", placeholder = "Filter by name…", `data-list` = "bf-list"),
              if (mode != "dir") checkboxInput(ns(id("all")), "Show all files", FALSE), checkboxInput(ns(id("hid")), "Show hidden", FALSE)),
          uiOutput(ns(id("list"))), uiOutput(ns(id("status"))))),                                   # (the modal is sent to the page right after this call returns)
      footer = tagList(actionButton(ns(id("cancel")), "Cancel", class = "btn-default"),
        if (mode %in% c("dir", "any")) actionButton(ns(id("pickdir")), "Use this folder", class = "btn-primary"),
        if (mode %in% c("file", "any")) actionButton(ns(id("pick")), if (multiple) "Choose the selected files" else "Choose this file", class = "btn-primary"))), session)
  }
  list(open = open)
}
