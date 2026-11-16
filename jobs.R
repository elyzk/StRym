# jobs.R -- background jobs for the GUI (a job keeps running if the browser tab is closed).
# State of each job lives in a small JSON file on the Linux disk, so a reopened page can reattach.

JOB_DIR <- path.expand("~/.tuxnet_r")
dir.create(JOB_DIR, showWarnings = FALSE)

# Is process `pid` alive (not a zombie) and really one of ours (command line contains `pattern`)?
pid_alive_pattern <- function(pid, pattern) {
  if (length(pid) != 1 || is.na(pid)) return(FALSE)
  f <- sprintf("/proc/%d/stat", as.integer(pid))
  if (dir.exists("/proc/self")) {                                       # Linux / WSL: read /proc directly (fast)
    if (!file.exists(f)) return(FALSE)
    if (sub(".*\\) (.).*", "\\1", readLines(f, warn = FALSE)[1]) == "Z") return(FALSE)
    cl <- tryCatch(readBin(sprintf("/proc/%d/cmdline", as.integer(pid)), "raw", 4096), error = function(e) raw())
    return(grepl(pattern, rawToChar(cl[cl != as.raw(0)]), fixed = TRUE))  # guards against pid reuse
  }
  out <- tryCatch(suppressWarnings(system2("ps", c("-p", as.integer(pid), "-o", "stat=,command="), stdout = TRUE, stderr = FALSE)),
                  error = function(e) character())                      # macOS (no /proc): ask ps
  if (!length(out) || !nzchar(trimws(out[1]))) return(FALSE)
  if (grepl("^\\s*Z", out[1])) return(FALSE)
  grepl(pattern, out[1], fixed = TRUE)
}

# Kill a process and everything it started (children first).
kill_tree <- function(pid) {
  ps <- tryCatch(system2("ps", c("-eo", "pid=,ppid="), stdout = TRUE), error = function(e) character())
  m <- do.call(rbind, lapply(strsplit(trimws(ps), "\\s+"), as.integer))
  kids <- function(p) { k <- m[m[, 2] == p, 1]; c(unlist(lapply(k, kids)), k) }
  for (p in c(if (!is.null(m)) kids(pid), pid)) tools::pskill(p, tools::SIGTERM)
}

job_state_file <- function(name) file.path(JOB_DIR, paste0(name, "_state.json"))
job_load <- function(name) {
  f <- job_state_file(name)
  if (file.exists(f)) tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
}
job_alive <- function(st, script) !is.null(st) && !is.null(st$pid) && pid_alive_pattern(st$pid, script)

# Start `Rscript <script> config.json` in the background. Returns the saved state.
job_start <- function(name, script, cfg, app_dir, out = NA_character_, result = NA_character_) {
  old <- job_load(name)                                                   # a second start (another browser tab, a double click) must not orphan a running job
  if (job_alive(old, script)) stop("A run of this kind is already in progress. Wait for it to finish or cancel it first.", call. = FALSE)
  cfgfile <- file.path(JOB_DIR, paste0(name, "_config.json"))
  logfile <- file.path(JOB_DIR, paste0(name, "_last_run.log"))
  jsonlite::write_json(cfg, cfgfile, auto_unbox = TRUE, na = "null", null = "null")
  p <- processx::process$new("Rscript", c(file.path(app_dir, script), cfgfile), stdout = logfile, stderr = "2>&1",
                             cleanup = FALSE, cleanup_tree = FALSE, env = Sys.getenv())
  st <- list(pid = p$get_pid(), log = logfile, out = out, result = result, script = script,
             started = format(Sys.time()))
  tmp <- paste0(job_state_file(name), ".tmp"); jsonlite::write_json(st, tmp, auto_unbox = TRUE, na = "null"); file.rename(tmp, job_state_file(name))   # written whole or not at all
  st
}

# "ok" | "failed" | "cancelled" from the log of a finished job (the CLI runners print DONE:/ERROR:).
job_outcome <- function(logfile) {
  txt <- if (file.exists(logfile)) paste(readLines(logfile, warn = FALSE), collapse = "\n") else ""
  if (grepl("\nDONE: |^DONE: ", txt)) "ok" else if (grepl("ERROR:", txt, fixed = TRUE)) "failed"
  else if (grepl("Execution halted|\nKilled|Segmentation fault|\nError", txt)) "failed"      # the process died on its own (a missing package, out of memory): not a cancel
  else "cancelled"
}

job_log_tail <- function(logfile, n = 400) {
  if (is.null(logfile) || is.na(logfile) || !file.exists(logfile)) return("")
  paste(tail(readLines(logfile, warn = FALSE), n), collapse = "\n")
}

# Last "[PROGRESS 40] label" marker in a log -> list(pct, label)
job_progress <- function(txt) {
  m <- regmatches(txt, gregexpr("\\[PROGRESS [0-9]+\\][^\n]*", txt))[[1]]
  if (!length(m)) return(list(pct = 0, label = "Starting…"))
  last <- tail(m, 1)
  list(pct = as.integer(sub("\\[PROGRESS ([0-9]+)\\].*", "\\1", last)),
       label = trimws(sub("^\\[PROGRESS [0-9]+\\]", "", last)))
}
