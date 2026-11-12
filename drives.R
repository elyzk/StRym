# drives.R -- find, mount and describe Windows drives from inside WSL.
#
# WSL only auto-mounts Windows drives that existed when WSL started. A hard drive plugged in
# later (typical for big sequencing data) is missing from /mnt until mounted with
#   mount -t drvfs D: /mnt/d
# These helpers let the GUI show every drive, mount the missing ones, and translate paths.

# Windows drives as seen by Windows (via interop), or NULL if interop is unavailable.
# Returns data.frame(letter, label, fs, size, free).
.ps_path <- function() {
  ps <- Sys.which("powershell.exe")
  if (!nzchar(ps)) ps <- "/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
  if (file.exists(ps)) ps else NULL
}
.VOLS_CMD <- paste0("Get-Volume | Where-Object { $_.DriveLetter } | ForEach-Object { ",
                    "'{0}|{1}|{2}|{3}|{4}' -f $_.DriveLetter, $_.FileSystemLabel, $_.FileSystem, ",
                    "$_.Size, $_.SizeRemaining }")
.windows_volumes <- function() {
  ps <- .ps_path(); if (is.null(ps)) return(NULL)
  out <- tryCatch(suppressWarnings(system2(ps, c("-NoProfile", "-NonInteractive", "-Command", shQuote(.VOLS_CMD)),
                                           stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  .parse_volumes(out)
}
.parse_volumes <- function(out) {
  out <- trimws(gsub("\r", "", out))
  out <- out[grepl("^[A-Za-z]\\|", out)]
  if (!length(out)) return(NULL)
  p <- strsplit(out, "|", fixed = TRUE)
  data.frame(letter = toupper(vapply(p, `[`, "", 1)), label = vapply(p, `[`, "", 2),
             fs = vapply(p, `[`, "", 3),
             size = suppressWarnings(as.numeric(vapply(p, `[`, "", 4))),
             free = suppressWarnings(as.numeric(vapply(p, `[`, "", 5))), stringsAsFactors = FALSE)
}

.is_mounted <- function(mp) {
  m <- tryCatch(readLines("/proc/mounts", warn = FALSE), error = function(e) character())
  any(vapply(strsplit(m, " ", fixed = TRUE), function(x) length(x) >= 2 && x[2] == mp, logical(1)))
}

# The /mnt/<letter> drive connections Ubuntu currently holds (whether or not Windows still has the drive).
.wsl_drive_mounts <- function() {
  m <- tryCatch(readLines("/proc/mounts", warn = FALSE), error = function(e) character())
  mp <- vapply(strsplit(m, " ", fixed = TRUE), function(x) if (length(x) >= 2) x[2] else "", "")
  sort(unique(mp[grepl("^/mnt/[a-z]$", mp)]))
}

# Does a mounted drive actually WORK? When a drive is unplugged or ejected and plugged back in, Ubuntu keeps the old connection ("mounted") but
# every read fails, so "mounted" alone is not enough. A quick `ls` with a time limit tells the difference. Answers are remembered for `ttl` seconds.
.usable_cache <- new.env()
.forget_usable <- function(mp = NULL) if (is.null(mp)) rm(list = ls(.usable_cache), envir = .usable_cache) else suppressWarnings(rm(list = mp, envir = .usable_cache))
.mount_usable <- function(mp, timeout = 4, ttl = 0) {
  now <- as.numeric(Sys.time()); hit <- .usable_cache[[mp]]
  if (ttl > 0 && !is.null(hit) && now - hit$t < ttl) return(hit$ok)
  r <- tryCatch(processx::run("ls", c("-A", mp), error_on_status = FALSE, timeout = timeout), error = function(e) NULL)
  ok <- !is.null(r) && identical(suppressWarnings(as.integer(r$status)), 0L)
  assign(mp, list(t = now, ok = ok), envir = .usable_cache); ok
}

# One row per Windows drive: letter, label, fs, size, free, mount point, mounted (logical), usable (logical) and state:
#   "ready"        connected and working
#   "not mounted"  Windows has the drive but Ubuntu is not connected to it yet
#   "stale"        connected, but the connection is dead (the drive was unplugged or ejected and plugged back in)
#   "removed"      Ubuntu is still connected to a drive that Windows no longer has
# Falls back to whatever is already mounted under /mnt if Windows can't be queried.
list_drives <- function() .assemble_drives(.windows_volumes(), check_usable = TRUE)
# check_usable = FALSE skips the `ls` test of every drive (up to 4 s each for a sleeping or unplugged one): "mounted" then counts as usable.
.assemble_drives <- function(vols, check_usable = TRUE) {
  windows_ok <- !is.null(vols)
  if (is.null(vols)) {
    cand <- setdiff(list.files("/mnt"), c("wsl", "wslg"))
    cand <- cand[nchar(cand) == 1]
    if (!length(cand)) return(data.frame(letter = character(), label = character(), fs = character(),
                                         size = numeric(), free = numeric(), mount = character(),
                                         mounted = logical(), usable = logical(), state = character(), stringsAsFactors = FALSE))
    vols <- data.frame(letter = toupper(cand), label = "", fs = "", size = NA_real_, free = NA_real_,
                       stringsAsFactors = FALSE)
  }
  vols$mount <- paste0("/mnt/", tolower(vols$letter))
  vols$mounted <- vapply(vols$mount, .is_mounted, logical(1))
  vols$usable <- vols$mounted & (if (check_usable) vapply(seq_len(nrow(vols)), function(i) vols$mounted[i] && .mount_usable(vols$mount[i]), logical(1)) else TRUE)
  vols$state <- ifelse(!vols$mounted, "not mounted", ifelse(vols$usable, "ready", "stale"))
  gone <- if (windows_ok) setdiff(.wsl_drive_mounts(), vols$mount) else character()
  if (length(gone)) vols <- rbind(vols, data.frame(letter = toupper(sub("^/mnt/", "", gone)), label = "", fs = "", size = NA_real_, free = NA_real_,
                                                   mount = gone, mounted = TRUE, usable = FALSE, state = "removed", stringsAsFactors = FALSE))
  vols[order(vols$letter), ]
}

# list_drives() remembered for a few seconds: the Browse… windows ask for the drive list at every click, and asking Windows takes a second or two.
.drives_cache <- new.env()
list_drives_cached <- function(ttl = 8) {
  now <- as.numeric(Sys.time()); hit <- .drives_cache$v
  if (!is.null(hit) && now - hit$t < ttl) return(hit$d)
  d <- list_drives(); .drives_cache$v <- list(t = now, d = d); d
}
refresh_drives <- function() { .drives_cache$v <- NULL; .roots_cache$v <- NULL; .forget_usable(); invisible(TRUE) }

# The drive list for the Browse windows, WITHOUT ever waiting. shinyFiles asks for the roots at every click in the window, and asking
# Windows (PowerShell) plus testing each drive takes seconds, during which the whole app is frozen and the page sits grey. So: the
# last answer is returned at once, and when it is older than `ttl` seconds Windows is asked again in a background process whose
# answer is picked up by a later call. The very first call (nothing known yet) lists what is already mounted under /mnt.
.roots_cache <- new.env(); .vols_proc <- new.env()
.refresh_drives_async <- function() {
  p <- .vols_proc$p; now <- as.numeric(Sys.time())
  if (!is.null(p) && p$is_alive()) { if (now - .vols_proc$t > 30) { try(p$kill(), silent = TRUE); .vols_proc$p <- NULL }; return(invisible()) }   # a hung PowerShell is stopped after 30 s
  if (!is.null(.vols_proc$t) && now - .vols_proc$t < 20) return(invisible())                      # asked a moment ago: never start one per call
  ps <- .ps_path(); if (is.null(ps)) return(invisible())
  .vols_proc$t <- now
  .vols_proc$p <- tryCatch(processx::process$new(ps, c("-NoProfile", "-NonInteractive", "-Command", .VOLS_CMD), stdout = "|", stderr = "|"),
                           error = function(e) NULL)
  invisible()
}
.collect_drives_async <- function() {
  p <- .vols_proc$p; if (is.null(p) || p$is_alive()) return(invisible())
  .vols_proc$p <- NULL
  out <- tryCatch(p$read_all_output_lines(), error = function(e) character())
  v <- .parse_volumes(out)
  if (!is.null(v)) .roots_cache$v <- list(t = as.numeric(Sys.time()), d = .assemble_drives(v, check_usable = FALSE))
  invisible()
}
drives_nowait <- function(ttl = 30) {
  .collect_drives_async()
  hit <- .roots_cache$v; now <- as.numeric(Sys.time())
  if (is.null(hit)) {
    full <- .drives_cache$v                                    # the blocking list, if something already asked for it
    .refresh_drives_async()
    return(if (!is.null(full)) full$d else .assemble_drives(NULL, check_usable = FALSE))
  }
  if (now - hit$t > ttl) .refresh_drives_async()
  hit$d
}

# Mount (or re-mount) a Windows drive letter. A dead connection (state "stale") is dropped first. Returns list(ok, message).
mount_drive <- function(letter) {
  letter <- toupper(trimws(letter))
  if (!grepl("^[A-Z]$", letter)) return(list(ok = FALSE, message = "Pick a drive letter."))
  mp <- paste0("/mnt/", tolower(letter)); refresh_drives()
  if (.is_mounted(mp)) {
    if (.mount_usable(mp)) return(list(ok = TRUE, message = paste0(letter, ": is already connected and working.")))
    suppressWarnings(system2("umount", c("-l", mp), stdout = TRUE, stderr = TRUE))      # the old connection is dead: let go of it, then mount fresh
    .forget_usable(mp)
  }
  dir.create(mp, showWarnings = FALSE)
  out <- suppressWarnings(system2("mount", c("-t", "drvfs", paste0(letter, ":"), mp), stdout = TRUE, stderr = TRUE))
  ok <- .is_mounted(mp) && .mount_usable(mp); refresh_drives()
  if (ok) return(list(ok = TRUE, message = paste0("Mounted ", letter, ": at ", mp, ". It is ready to use: type a path such as ", letter, ":\\data in any box, or use Browse…")))
  why <- trimws(paste(out, collapse = " "))
  list(ok = FALSE, message = paste0(
    if (grepl("does not exist", why, fixed = TRUE)) sprintf("Windows has no drive %s: right now. Check that the drive is plugged in and has a drive letter (Windows Disk Management), then press Refresh. ", letter)
    else if (Sys.info()[["user"]] != "root") sprintf("Could not mount %s: (StRym must run as root inside Ubuntu, or run  sudo mount -t drvfs %s: %s ). ", letter, letter, mp)
    else sprintf("Could not mount %s: %s ", letter, if (nzchar(why)) paste0("(", substr(why, 1, 160), ")") else ""),
    "Close any File Explorer window showing the drive, press Refresh and try Mount drive again; if it still fails, see “If it still does not work” below."))
}

# Let go of a drive so Windows can eject it. The caller makes sure no StRym run is using it. Returns list(ok, message).
unmount_drive <- function(letter) {
  letter <- toupper(trimws(letter))
  if (!grepl("^[A-Z]$", letter)) return(list(ok = FALSE, message = "Pick a drive letter."))
  if (letter == "C") return(list(ok = FALSE, message = "C: is your Windows system drive: StRym keeps it connected."))
  mp <- paste0("/mnt/", tolower(letter))
  if (!.is_mounted(mp)) { refresh_drives(); return(list(ok = TRUE, message = paste0(letter, ": is not connected to StRym. You can eject it in Windows."))) }
  out <- suppressWarnings(system2("umount", mp, stdout = TRUE, stderr = TRUE))
  if (.is_mounted(mp)) out <- c(out, suppressWarnings(system2("umount", c("-l", mp), stdout = TRUE, stderr = TRUE)))   # busy: detach it now
  refresh_drives(); ok <- !.is_mounted(mp)
  list(ok = ok, message = if (ok) paste0("Released ", letter, ":. You can eject it in Windows now. To use it in StRym again, plug it in and press Mount drive.")
                          else paste0("Could not release ", letter, ": ", trimws(paste(out, collapse = " "))))
}

# Named vector of roots for the shinyFiles pickers. Every Windows drive that is there is included, mounted or not: the
# /mnt/<letter> folder is created up front so the picker works as soon as the drive is mounted. (The app hands the pickers a FUNCTION that
# calls this, so a drive mounted later shows up in the next Browse… window without reloading the page.)
.root_names <- new.env()          # a drive keeps the name it was first shown under: the window refers to a root by its name, so it must not change
drive_roots <- function(drives = drives_nowait()) {
  r <- c(Home = path.expand("~"))
  if ("state" %in% names(drives)) drives <- drives[drives$state != "removed", , drop = FALSE]
  if (nrow(drives)) {
    for (mp in drives$mount) dir.create(mp, showWarnings = FALSE)
    nm <- ifelse(nzchar(drives$label), sprintf("%s: %s", drives$letter, drives$label),
                 paste0(drives$letter, ":"))
    for (i in seq_along(nm)) { mp <- drives$mount[i]; if (is.null(.root_names[[mp]])) .root_names[[mp]] <- nm[i]; nm[i] <- .root_names[[mp]] }
    r <- c(r, setNames(drives$mount, nm))
  }
  c(r, Root = "/")
}

# "/mnt/d/foo/bar" -> "D:\foo\bar" (for display); other paths unchanged.
to_windows_path <- function(p) {
  m <- regmatches(p, regexec("^/mnt/([a-z])(/.*)?$", p))[[1]]
  if (length(m)) paste0(toupper(m[2]), ":", gsub("/", "\\\\", if (is.na(m[3]) || m[3] == "") "\\" else m[3])) else p
}

# Translate a user-typed path: "D:\data\x" -> "/mnt/d/data/x", strip quotes/whitespace.
clean_path <- function(p) {
  p <- gsub("^[\"']|[\"']$", "", trimws(p))
  p <- sub("^file:///", "", p)                                                                   # a pasted file:// address
  unc <- "^\\\\\\\\wsl(\\.localhost|\\$)\\\\[^\\\\]+"                                              # \\wsl$\Ubuntu\home\me (Explorer's Linux view) -> /home/me
  if (grepl(unc, p)) { p <- gsub("\\\\", "/", sub(unc, "", p)); if (!nzchar(p)) p <- "/" }
  if (grepl("^[A-Za-z]:([\\\\/]|$)", p))
    p <- paste0("/mnt/", tolower(substr(p, 1, 1)), "/", gsub("\\\\", "/", substring(p, 4)))
  if (identical(p, "/")) return(p)                                                               # the root folder is not an empty path
  sub("/+$", "", p)
}

# If `p` is under /mnt/<letter> but that drive isn't mounted, return a helpful message.
unmounted_drive_problem <- function(p) {
  m <- regmatches(p, regexec("^/mnt/([a-z])(/|$)", p))[[1]]
  if (!length(m)) return("")
  mp <- paste0("/mnt/", m[2])
  if (!.is_mounted(mp))
    return(sprintf("Drive %s: is not connected to StRym yet. Open the Drives tab, press Refresh, then “Mount drive”.", toupper(m[2])))
  if (!.mount_usable(mp, ttl = 15))
    return(sprintf("Drive %s: is connected but not responding (it was probably unplugged or ejected and plugged back in). Open the Drives tab and press “Mount drive” to reconnect it.", toupper(m[2])))
  ""
}

# ---- finding files that were dragged and dropped ---------------------------------------------------------------
# A browser only tells us the NAME (and size) of a dropped file, never where it lives. So StRym looks for it: first in the
# folders the researcher is already working in (`hints`), then in the usual places (Downloads, Desktop, Documents, the
# user profile, every mounted drive), shallow folders first, within a time budget. A file only counts as found when its
# size matches too. Returns a named character vector (names = `names`): the path, or NA when not found.
locate_dropped <- function(names, sizes = rep(NA_real_, length(names)), dirs = rep(FALSE, length(names)),
                           hints = character(), max_seconds = 45) {
  res <- setNames(rep(NA_character_, length(names)), names)
  sizes <- suppressWarnings(as.numeric(sizes)); dirs <- as.logical(dirs)
  ok <- function(p, i) if (isTRUE(dirs[i])) dir.exists(p)
    else file.exists(p) && !dir.exists(p) && (is.na(sizes[i]) || isTRUE(file.size(p) == sizes[i]))
  hints <- unique(hints[nzchar(hints)]); hints <- hints[dir.exists(hints)]
  hints <- unique(c(hints, dirname(hints)))
  for (i in seq_along(names)) for (h in hints) { p <- file.path(h, names[i]); if (ok(p, i)) { res[i] <- p; break } }
  if (!anyNA(res)) return(res)

  users <- list.dirs("/mnt/c/Users", recursive = FALSE, full.names = TRUE)
  users <- users[!basename(users) %in% c("Public", "Default", "Default User", "All Users")]
  drives <- setdiff(list.dirs("/mnt", recursive = FALSE, full.names = TRUE), c("/mnt/c", "/mnt/wsl", "/mnt/wslg"))
  drives <- drives[nchar(basename(drives)) == 1]
  places <- list()                                                         # list(dir, depth), most likely first
  for (u in users) for (f in c("Downloads", "Desktop", "Documents")) places[[length(places) + 1]] <- list(file.path(u, f), 6)
  for (u in users) places[[length(places) + 1]] <- list(u, 4)
  for (dv in drives) places[[length(places) + 1]] <- list(dv, 4)
  places[[length(places) + 1]] <- list(path.expand("~"), 4)
  places[[length(places) + 1]] <- list("/mnt/c", 3)
  t0 <- Sys.time()
  for (pl in places) {
    todo <- which(is.na(res)); if (!length(todo)) break
    left <- max_seconds - as.numeric(difftime(Sys.time(), t0, units = "secs")); if (left <= 1) break
    if (!dir.exists(pl[[1]])) next
    esc <- gsub("([][*?\\])", "\\\1", names[todo])
    pat <- unlist(lapply(esc, function(n) c("-o", "-name", n)))[-1]
    prune <- c("(", "-name", "AppData", "-o", "-name", "node_modules", "-o", "-name", "$RECYCLE.BIN", "-o", "-name", "System Volume Information", ")", "-prune", "-o")
    out <- tryCatch(processx::run("find", c(pl[[1]], "-maxdepth", pl[[2]], prune, "(", pat, ")", "-print"),
                                  error_on_status = FALSE, timeout = min(20, left))$stdout, error = function(e) "")
    cand <- strsplit(out, "\n", fixed = TRUE)[[1]]; cand <- cand[nzchar(cand)]
    for (i in todo) {
      hit <- cand[basename(cand) == names[i]]; hit <- hit[vapply(hit, ok, logical(1), i = i)]
      if (length(hit)) res[i] <- hit[order(nchar(hit))][1]
    }
  }
  res
}
