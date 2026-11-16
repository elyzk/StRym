# drives.R -- find and describe the drives (external disks, USB sticks) StRym can read from and write to, on every platform.
#
# Windows (StRym runs inside WSL): WSL only connects the Windows drives that existed when it started. A hard drive plugged in later
#   (typical for big sequencing data) is missing from /mnt until it is mounted with  mount -t drvfs D: /mnt/d . On this platform the
#   helpers below list every Windows drive, mount the missing ones, release them before ejecting, and translate D:\path <-> /mnt/d/path.
# macOS and Linux: the operating system mounts drives by itself (/Volumes/<name> on a Mac, /media/<user>/<name> or /run/media/<user>/<name>
#   on most Linux systems, /mnt/<name> on others), so StRym only lists what is mounted and checks that it answers.
# strym_platform() tells the two worlds apart; everything that is specific to WSL is skipped elsewhere.

# "wsl" (Linux inside Windows), "linux", "mac" or "windows" (R running natively on Windows). options(strym.platform = ...) or the environment variable STRYM_PLATFORM override it (tests, or to look at another platform's Drives page).
strym_platform <- function() {
  o <- getOption("strym.platform"); if (!is.null(o)) return(as.character(o)[1])
  e <- Sys.getenv("STRYM_PLATFORM"); if (e %in% c("wsl", "linux", "mac", "windows")) return(e)                    # for looking at the other platforms' pages (and for tests)
  if (identical(.Platform$OS.type, "windows")) return("windows")
  if (identical(Sys.info()[["sysname"]], "Darwin")) return("mac")
  v <- tryCatch(readLines("/proc/version", n = 1, warn = FALSE), error = function(e) character())
  if (length(v) && grepl("microsoft|wsl", v[1], ignore.case = TRUE)) "wsl" else "linux"
}
strym_is_wsl <- function() identical(strym_platform(), "wsl")

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

# ---- macOS, Linux and native Windows: the operating system mounts the drives; StRym lists them ---------------------------------------------
# Free and total size of the file system that holds `path` (bytes), NA when df cannot tell.
.df_space <- function(path) {
  r <- tryCatch(processx::run("df", c("-Pk", path), error_on_status = FALSE, timeout = 5)$stdout, error = function(e) "")
  l <- strsplit(r, "\n", fixed = TRUE)[[1]]; l <- l[nzchar(l)]
  if (length(l) < 2) return(c(size = NA_real_, free = NA_real_))
  f <- strsplit(trimws(l[length(l)]), "[[:space:]]+")[[1]]
  c(size = suppressWarnings(as.numeric(f[length(f) - 4])) * 1024, free = suppressWarnings(as.numeric(f[length(f) - 2])) * 1024)
}
# The mount points that are external or extra drives. Mac: /Volumes/<name> (the start disk, a link to /, is left out). Linux: the usual places of
# the desktop auto-mounters and /mnt, for real devices only (no proc, tmpfs, snap loops ...). Native Windows: the drive letters that exist.
.volume_mounts <- function(platform = strym_platform(), proc_mounts = NULL) {          # proc_mounts: the lines of /proc/mounts (a test passes its own)
  if (identical(platform, "mac")) {
    d <- list.dirs("/Volumes", recursive = FALSE, full.names = TRUE)
    return(d[!nzchar(Sys.readlink(d))])
  }
  if (identical(platform, "windows")) { l <- paste0(LETTERS, ":/"); return(l[dir.exists(l)]) }
  m <- if (!is.null(proc_mounts)) proc_mounts else tryCatch(readLines("/proc/mounts", warn = FALSE), error = function(e) character())
  f <- strsplit(m, " ", fixed = TRUE); f <- f[lengths(f) >= 3]
  dev <- vapply(f, `[`, "", 1); mp <- vapply(f, `[`, "", 2); fs <- vapply(f, `[`, "", 3)
  mp <- gsub("\\040", " ", mp, fixed = TRUE)                                              # /proc/mounts writes a space as \040
  ok <- grepl("^(/media/|/run/media/|/mnt/|/Volumes/)", mp) & (grepl("^/dev/", dev) | fs %in% c("fuseblk", "exfat", "ntfs", "ntfs3", "vfat", "cifs", "nfs", "nfs4", "9p")) &
        !grepl("^/dev/loop", dev)
  sort(unique(mp[ok]))
}
# Same columns as .assemble_drives() so every consumer works on either: `letter` holds the drive name (the folder name), `state` is ready / stale.
.generic_drives <- function(platform = strym_platform(), check_usable = TRUE) {
  mp <- .volume_mounts(platform)
  if (!length(mp)) return(data.frame(letter = character(), label = character(), fs = character(), size = numeric(), free = numeric(), mount = character(),
                                     mounted = logical(), usable = logical(), state = character(), stringsAsFactors = FALSE))
  nm <- sub(":$", "", basename(sub("/+$", "", mp)))
  sp <- vapply(mp, .df_space, c(size = 0, free = 0)); sp <- matrix(sp, nrow = 2, dimnames = list(c("size", "free"), NULL))
  us <- if (check_usable) vapply(mp, function(x) .mount_usable(x), logical(1)) else rep(TRUE, length(mp))
  data.frame(letter = nm, label = nm, fs = "", size = unname(sp["size", ]), free = unname(sp["free", ]), mount = mp, mounted = TRUE, usable = unname(us),
             state = ifelse(us, "ready", "stale"), stringsAsFactors = FALSE)
}

# One row per Windows drive: letter, label, fs, size, free, mount point, mounted (logical), usable (logical) and state:
#   "ready"        connected and working
#   "not mounted"  Windows has the drive but Ubuntu is not connected to it yet
#   "stale"        connected, but the connection is dead (the drive was unplugged or ejected and plugged back in)
#   "removed"      Ubuntu is still connected to a drive that Windows no longer has
# Falls back to whatever is already mounted under /mnt if Windows can't be queried.
list_drives <- function() if (strym_is_wsl()) .assemble_drives(.windows_volumes(), check_usable = TRUE) else .generic_drives(check_usable = TRUE)
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
  vols$usable <- vols$mounted & (if (check_usable) vapply(seq_len(nrow(vols)), function(i) vols$mounted[i] && .mount_usable(vols$mount[i], timeout = 2), logical(1)) else TRUE)
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
refresh_drives <- function() {
  .drives_cache$v <- NULL; .roots_cache$v <- NULL; .forget_usable()
  if (strym_is_wsl()) { p <- .vols_proc$p; if (is.null(p) || !p$is_alive()) { .vols_proc$t <- 0; .refresh_drives_async() } }       # ask Windows again IN THE BACKGROUND (unless the question is already on its way); the old answer stays on screen until the new one arrives
  invisible(TRUE)
}

# The drive list for Browse, returned without ever waiting. Asking Windows through PowerShell and testing every drive can
# take seconds, so callers get the last answer immediately. When it is stale, a background process refreshes it for a
# later call. The first call, before Windows has answered, lists mounts already visible under /mnt.
.roots_cache <- new.env(); .vols_proc <- new.env(); .vols_last <- new.env()          # .vols_last$v: the last answer of Windows (data.frame of volumes)
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
  if (!is.null(v)) { .vols_last$v <- v; .roots_cache$v <- list(t = as.numeric(Sys.time()), d = .assemble_drives(v, check_usable = FALSE)) }
  invisible()
}
drives_nowait <- function(ttl = 30) {
  if (!strym_is_wsl()) return(.generic_drives(check_usable = FALSE))                      # nothing to ask Windows: reading /proc/mounts or /Volumes is instant
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

# The list for the Drives page. It never waits for Windows: it is built at once from the last answer Windows gave (or from what is mounted under /mnt, before
# the first answer), and the page asks again until the background question has finished (drives_pending()). Unlike drives_nowait() it also tests whether each
# mounted drive still answers (a quick `ls`, at most 2 s per drive), so "needs reconnecting" is shown.
drives_for_table <- function() {
  if (!strym_is_wsl()) return(.generic_drives(check_usable = TRUE))
  .collect_drives_async()
  if (is.null(.vols_proc$p) && is.null(.vols_last$v)) .refresh_drives_async()
  .assemble_drives(.vols_last$v, check_usable = TRUE)
}
drives_pending <- function() strym_is_wsl() && ((!is.null(.vols_proc$p) && .vols_proc$p$is_alive()) || is.null(.vols_last$v))

# Mount (or re-mount) a Windows drive letter. A dead connection (state "stale") is dropped first. Returns list(ok, message).
mount_drive <- function(letter) {
  if (!strym_is_wsl()) return(list(ok = FALSE, message = "Your operating system mounts drives by itself: plug the drive in and press Refresh."))
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

# Which programs still use a mounted drive? A program uses it when its working folder is on the drive or it holds a file open there. This is what
# makes "target is busy" (and Windows' "this device is in use"): the list names them so the cause is visible instead of guessed.
# Returns data.frame(pid, command, strym) -- strym = TRUE for programs StRym started itself (analysis runs, the live browser, the tools they call).
drive_holders <- function(mp) {
  mp <- sub("/+$", "", mp); me <- Sys.getpid(); pids <- list.files("/proc", "^[0-9]+$"); pids <- pids[as.integer(pids) != me]
  rows <- lapply(pids, function(p) {
    cwd <- tryCatch(Sys.readlink(file.path("/proc", p, "cwd")), error = function(e) "")
    fd <- tryCatch(list.files(file.path("/proc", p, "fd"), full.names = TRUE), error = function(e) character()); fl <- if (length(fd)) Sys.readlink(fd) else character()
    all <- c(cwd, fl); if (!any(startsWith(paste0(all, "/"), paste0(mp, "/")))) return(NULL)
    raw <- tryCatch(readBin(file.path("/proc", p, "cmdline"), "raw", 3000), error = function(e) raw())
    cmd <- if (length(raw)) trimws(rawToChar(replace(raw, raw == as.raw(0), as.raw(32)))) else ""                   # the arguments are separated by NUL bytes: show them with spaces
    data.frame(pid = as.integer(p), command = substr(cmd, 1, 110), stringsAsFactors = FALSE)
  })
  h <- do.call(rbind, rows); if (is.null(h)) return(data.frame(pid = integer(), command = character(), strym = logical(), stringsAsFactors = FALSE))
  h$strym <- grepl("run_(analysis|genist|rtpstar)_cli\\.R|run_cli\\.R|runApp|cuffdiff|cufflinks|cuffmerge|hisat2|samtools|minimap2|featureCounts|fastq-mcf|pheatmap|Rscript", h$command)
  h
}

# Let go of a drive so Windows can eject it (Windows / WSL only; elsewhere the system ejects). The caller makes sure no StRym run is going. Steps:
#   1. write everything out (sync);  2. stop StRym's own helpers that still sit on the drive (a live browser started from the project folder keeps it busy);
#   3. unmount; when it is busy, name the programs that hold it and, if they are not StRym's, detach it anyway ("lazy"), saying so;
#   4. ask Windows to flush its own write cache for the drive (best effort).
# Returns list(ok, message, holders).
unmount_drive <- function(letter) {
  if (!strym_is_wsl()) return(list(ok = FALSE, message = "Eject the drive with your operating system (Finder, the file manager or Safely Remove Hardware) once no StRym run is using it.", holders = NULL))
  letter <- toupper(trimws(letter)); res <- function(ok, msg, h = NULL) list(ok = ok, message = msg, holders = h)
  if (!grepl("^[A-Z]$", letter)) return(res(FALSE, "Pick a drive letter."))
  if (letter == "C") return(res(FALSE, "C: is your Windows system drive: StRym keeps it connected."))
  mp <- paste0("/mnt/", tolower(letter))
  if (!.is_mounted(mp)) { refresh_drives(); return(res(TRUE, paste0(letter, ": is not connected to StRym. You can eject it in Windows."))) }
  suppressWarnings(system2("sync", stdout = FALSE, stderr = FALSE))
  h <- drive_holders(mp)
  jobs <- h[grepl("run_(analysis|genist|rtpstar)_cli\\.R|run_cli\\.R", h$command), , drop = FALSE]                  # a StRym run that is still going must finish (or be cancelled) first
  if (nrow(jobs)) return(res(FALSE, paste0("A StRym run is still using ", letter, ": (", sub("^.*/", "", sub(" .*$", "", jobs$command[1])), " ", sub("^\\S+\\s+", "", jobs$command[1]), "). Wait for it to finish, or cancel it, then release the drive."), h))
  live <- h[grepl("runApp", h$command), , drop = FALSE]                                      # the live results browser: StRym's own, safe to stop
  for (pid in live$pid) tools::pskill(pid, tools::SIGTERM)
  if (nrow(live)) { Sys.sleep(1); h <- drive_holders(mp) }
  out <- suppressWarnings(system2("umount", mp, stdout = TRUE, stderr = TRUE))
  if (.is_mounted(mp)) { Sys.sleep(2); out <- suppressWarnings(system2("umount", mp, stdout = TRUE, stderr = TRUE)); h <- drive_holders(mp) }  # a program that was just finishing
  lazy <- FALSE
  if (.is_mounted(mp) && nrow(h)) {                                                           # still busy and we can name who: do NOT pretend it is free. Leave it connected so that pressing the button again, after they are closed, really releases it
    names <- unique(sub("^.*/", "", sub(" .*$", "", h$command)))
    return(res(FALSE, paste0("Not ready to eject yet: ", paste(names, collapse = ", "), if (nrow(h) == 1) " is" else " are", " still using ", letter, ":. Close ", if (nrow(h) == 1) "it" else "them",
                             " (an Ubuntu window whose folder is on the drive is a common cause), then press Release drive again."), h))
  }
  if (.is_mounted(mp)) { lazy <- TRUE; out <- c(out, suppressWarnings(system2("umount", c("-l", mp), stdout = TRUE, stderr = TRUE))) }               # busy for a reason we cannot see: detach it now, the system lets go when the last user does
  refresh_drives(); ok <- !.is_mounted(mp)
  if (!ok) return(res(FALSE, paste0("Could not release ", letter, ": ", trimws(paste(out, collapse = " "))), h))
  ps <- .ps_path()                                                                           # Windows keeps its own write cache: flush it so the drive can be removed safely (needs no extra rights on most systems)
  if (!is.null(ps)) try(processx::run(ps, c("-NoProfile", "-NonInteractive", "-Command", sprintf("Write-VolumeCache -DriveLetter %s", letter)), error_on_status = FALSE, timeout = 15), silent = TRUE)
  msg <- if (lazy) paste0("Released ", letter, ": (detached; something is still closing it). Wait a few seconds before ejecting.")
         else paste0("Released ", letter, ":. You can eject it in Windows now.")
  res(TRUE, paste0(msg, " If Windows still says the drive is in use, close File Explorer windows that show it and wait a little (the antivirus may be scanning new files)."), h)
}

# "/mnt/d/foo/bar" -> "D:\foo\bar" (for display); other paths unchanged.
to_windows_path <- function(p) {
  if (!strym_is_wsl()) return(p)
  m <- regmatches(p, regexec("^/mnt/([a-z])(/.*)?$", p))[[1]]
  if (length(m)) paste0(toupper(m[2]), ":", gsub("/", "\\\\", if (is.na(m[3]) || m[3] == "") "\\" else m[3])) else p
}

# Translate a user-typed path: "D:\data\x" -> "/mnt/d/data/x", strip quotes/whitespace.
clean_path <- function(p) {
  p <- gsub("^[\"']|[\"']$", "", trimws(p))
  p <- sub("^file:///", "", p)                                                                   # a pasted file:// address
  unc <- "^\\\\\\\\wsl(\\.localhost|\\$)\\\\[^\\\\]+"                                              # \\wsl$\Ubuntu\home\me (Explorer's Linux view) -> /home/me
  if (strym_is_wsl() && grepl(unc, p)) { p <- gsub("\\\\", "/", sub(unc, "", p)); if (!nzchar(p)) p <- "/" }
  if (strym_is_wsl() && grepl("^[A-Za-z]:([\\\\/]|$)", p))
    p <- paste0("/mnt/", tolower(substr(p, 1, 1)), "/", gsub("\\\\", "/", substring(p, 4)))
  if (identical(p, "/")) return(p)                                                               # the root folder is not an empty path
  sub("/+$", "", p)
}

# If `p` is under /mnt/<letter> but that drive isn't mounted, return a helpful message.
unmounted_drive_problem <- function(p) {
  if (!strym_is_wsl()) {                                                                   # Mac / Linux: a path under a drive that no longer answers
    mp <- .volume_mounts(); mp <- mp[nzchar(mp)]; if (!length(mp)) return("")                # (paste0(character(0), "/") is "/": an empty list must not match everything)
    hit <- mp[startsWith(paste0(sub("/+$", "", p), "/"), paste0(sub("/+$", "", mp), "/"))]
    if (length(hit)) { h <- hit[which.max(nchar(hit))]
      if (!.mount_usable(h, ttl = 15)) return(sprintf("The drive %s is not responding (it was probably unplugged or ejected and plugged back in). Plug it in again, open the Drives tab and press Refresh.", basename(h))) }
    return("")
  }
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

  wsl <- strym_is_wsl()
  user_root <- if (wsl) "/mnt/c/Users" else if (identical(strym_platform(), "mac")) "/Users" else "/home"
  users <- if (dir.exists(user_root)) list.dirs(user_root, recursive = FALSE, full.names = TRUE) else character()
  users <- users[!basename(users) %in% c("Public", "Default", "Default User", "All Users", "Shared", "Guest")]
  drives <- if (wsl) { d <- setdiff(list.dirs("/mnt", recursive = FALSE, full.names = TRUE), c("/mnt/c", "/mnt/wsl", "/mnt/wslg")); d[nchar(basename(d)) == 1] } else .volume_mounts()
  places <- list()                                                         # list(dir, depth), most likely first
  for (u in users) for (f in c("Downloads", "Desktop", "Documents")) places[[length(places) + 1]] <- list(file.path(u, f), 6)
  for (u in users) places[[length(places) + 1]] <- list(u, 4)
  for (dv in drives) places[[length(places) + 1]] <- list(dv, 4)
  places[[length(places) + 1]] <- list(path.expand("~"), 4)
  if (wsl) places[[length(places) + 1]] <- list("/mnt/c", 3)
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
