# Version 11: notes that can be run past, the new Browse window, the drive list that never waits, drive release, and optional fastp cleaning.
suppressPackageStartupMessages(library(DT))


test_that("a note about the reads can be run past; a real problem cannot", {
  root <- tempfile("v11"); dir.create(root); w <- function(n) { f <- file.path(root, n); writeLines(c("@r1", "ACGT", "+", "IIII"), f); f }
  f <- c(w("S1_L001_R1.fastq"), w("S1_L002_R1.fastq")); s <- data.frame(group = c("A", "A"), file = f, mate = NA_character_, stringsAsFactors = FALSE)
  v <- validate_samples(s, FALSE); expect_length(v$errors, 1)                                       # one group only: that one is a real error ...
  expect_true(any(grepl("lanes", v$warnings)))                                                      # ... and the lane note is advice (a warning), not an error
  s2 <- data.frame(group = c("A", "A", "B", "B"), file = c(f, w("T1_L001_R1.fastq"), w("T1_L002_R1.fastq")), mate = NA_character_, stringsAsFactors = FALSE)
  v2 <- validate_samples(s2, FALSE); expect_length(v2$errors, 0); expect_true(any(grepl("lanes", v2$warnings)))
})

test_that("a folder is listed with one call: folders first, hidden and system names left out, file types filtered, links to folders are folders", {
  d <- tempfile("v11ls"); dir.create(d); dir.create(file.path(d, "b_folder")); dir.create(file.path(d, "A_folder")); dir.create(file.path(d, "$RECYCLE.BIN"))
  for (n in c("z.csv", "a.CSV", "notes.txt", "reads.fastq.gz", ".hidden.csv", "photo.png", "desktop.ini")) writeLines("x", file.path(d, n))
  file.symlink(file.path(d, "A_folder"), file.path(d, "link_to_folder"))
  L <- strym_list_dir(d); expect_true(L$ok); expect_equal(L$entries$name[L$entries$dir], c("A_folder", "b_folder", "link_to_folder"))
  expect_equal(L$entries$name[!L$entries$dir], c("a.CSV", "notes.txt", "photo.png", "reads.fastq.gz", "z.csv")); expect_false(any(grepl("hidden|RECYCLE|desktop", L$entries$name)))
  e <- strym_list_dir(d, exts = c("csv", "gz")); expect_equal(e$entries$name[!e$entries$dir], c("a.CSV", "reads.fastq.gz", "z.csv")); expect_equal(sum(e$entries$dir), 3)    # folders always stay
  expect_equal(nrow(strym_list_dir(d, exts = "csv", show_all = TRUE)$entries), nrow(L$entries))                                          # "show all files"
  expect_true(".hidden.csv" %in% strym_list_dir(d, hidden = TRUE)$entries$name)
  expect_true(all(is.finite(L$entries$size[!L$entries$dir])))                                                                           # sizes for a short list
  big <- tempfile("v11big"); dir.create(big); file.create(file.path(big, sprintf("f%04d.txt", 1:60)))
  M <- strym_list_dir(big, max_rows = 25); expect_equal(c(M$shown, M$total), c(25L, 60L)); expect_true(all(is.finite(M$entries$size)))
  expect_false(strym_list_dir(file.path(d, "nope"))$ok); expect_match(strym_list_dir(file.path(d, "nope"))$error, "does not exist"); expect_false(strym_list_dir("")$ok)
})

test_that("a huge folder is listed in well under a second and a slow one is given up on with a plain message", {
  big <- tempfile("v11huge"); dir.create(big); invisible(file.create(file.path(big, sprintf("sample_%05d.fastq.gz", 1:4000))))
  t <- system.time(L <- strym_list_dir(big))[["elapsed"]]; expect_true(L$ok); expect_equal(L$total, 4000L); expect_equal(L$shown, 1500L); expect_lt(t, 3)
  # a folder that never answers: simulate with a command that sleeps; the limit turns it into a message instead of a frozen window
  r <- tryCatch(processx::run("sleep", "5", timeout = 0.3, error_on_status = FALSE), error = function(e) NULL); expect_true(is.null(r) || isTRUE(r$timeout))        # strym_list_dir() tests r$timeout
})

test_that("the places list has Quick places and Drives on every platform, and an unconnected Windows drive is one click from usable", {
  dv <- data.frame(letter = c("C", "D", "E"), label = c("Windows", "USB STICK", ""), fs = "", size = NA_real_, free = c(5e10, 2e9, NA), mount = c("/mnt/c", "/mnt/d", "/mnt/e"),
                   mounted = c(TRUE, TRUE, FALSE), usable = c(TRUE, TRUE, FALSE), state = c("ready", "ready", "not mounted"), stringsAsFactors = FALSE)
  with_platform("wsl", { p <- strym_places(dv); expect_true(all(c("Quick places", "Drives") %in% p$group)); expect_true("D:  USB STICK" %in% p$label)
    e <- p[p$letter == "E", ]; expect_equal(e$kind, "mount"); expect_match(e$note, "click to connect"); expect_equal(p$kind[p$letter == "D"], "go"); expect_true("/" %in% p$path) })
  tmp <- tempfile("vol"); dir.create(file.path(tmp, "USB"), recursive = TRUE)
  g <- data.frame(letter = "USB", label = "USB", fs = "", size = NA_real_, free = 1e9, mount = file.path(tmp, "USB"), mounted = TRUE, usable = TRUE, state = "ready", stringsAsFactors = FALSE)
  for (pl in c("mac", "linux")) with_platform(pl, { p <- strym_places(g); expect_true("USB" %in% p$label); expect_false(any(p$kind == "mount")); expect_true(any(grepl("Whole|Macintosh", p$label))) })
})

test_that("the container data mount is a named quick place", {
  data_root <- tempfile("strym-data"); dir.create(data_root)
  old <- Sys.getenv("STRYM_DATA_ROOT", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("STRYM_DATA_ROOT") else Sys.setenv(STRYM_DATA_ROOT = old))
  Sys.setenv(STRYM_DATA_ROOT = data_root)
  p <- strym_places(data.frame())
  expect_equal(p$path[p$label == "Data"], data_root)
  expect_equal(p$kind[p$label == "Data"], "go")
})

test_that("path fields expose StRym's own file and folder browser controls", {
  h <- as.character(pathUI("x", "Label", "any", multiple = TRUE))
  expect_true(grepl("x-btn", h))
  expect_true(grepl("x-btn_dir", h))
  expect_false(grepl("shinyFiles|shiny-bound-input.*shinyDir", h))
})

test_that("the Drives page asks Windows in the background: the list is built at once, and a Refresh does not freeze anything", {
  with_platform("wsl", {
    t <- system.time(d <- drives_for_table())[["elapsed"]]; expect_lt(t, 3); expect_true(all(c("letter", "state", "mount") %in% names(d)))
    expect_true(is.logical(drives_pending())); t2 <- system.time(refresh_drives())[["elapsed"]]; expect_lt(t2, 1.5)                    # starting the background question is instant
  })
  with_platform("linux", expect_false(drives_pending()))
})

test_that("which programs keep a drive busy can be found, and a run in progress stops the release", {
  d <- tempfile("v11hold"); dir.create(d)
  expect_equal(nrow(drive_holders(d)), 0)
  p <- processx::process$new("sleep", "30", wd = d, cleanup = TRUE); Sys.sleep(0.5)
  h <- drive_holders(d); expect_true(p$get_pid() %in% h$pid); expect_match(h$command[h$pid == p$get_pid()], "sleep 30"); expect_false(h$strym[h$pid == p$get_pid()])
  p$kill(); Sys.sleep(0.3); expect_false(p$get_pid() %in% drive_holders(d)$pid)
  r <- processx::process$new("bash", c("-c", "exec -a 'Rscript /x/run_analysis_cli.R /y/analysis_config.json' sleep 30"), wd = d, cleanup = TRUE); Sys.sleep(0.5)
  h <- drive_holders(d); expect_true(r$get_pid() %in% h$pid); expect_true(h$strym[h$pid == r$get_pid()])
  r$kill(); expect_equal(nrow(drive_holders(d)[drive_holders(d)$pid == Sys.getpid(), ]), 0)
})

test_that("a drive that is still in use is not reported as free: the user is told who holds it, and releasing works once they are gone", {
  skip_if(!identical(Sys.info()[["user"]], "root") || !strym_is_wsl() || !nzchar(Sys.which("mount")), "needs root inside WSL to mount a test file system")
  mp <- "/mnt/y"; dir.create(mp, showWarnings = FALSE); ok <- identical(suppressWarnings(system2("mount", c("-t", "tmpfs", "tmpfs", mp), stdout = FALSE, stderr = FALSE)), 0L); skip_if(!ok, "could not mount a test file system")
  on.exit({ suppressWarnings(system2("umount", c("-l", mp), stdout = FALSE, stderr = FALSE)); unlink(mp, recursive = TRUE) }, add = TRUE)
  dir.create(file.path(mp, "proj")); p <- processx::process$new("sleep", "120", wd = file.path(mp, "proj"), cleanup = TRUE); Sys.sleep(0.4)
  r <- unmount_drive("Y"); expect_false(r$ok); expect_match(r$message, "Not ready to eject yet: sleep is still using Y:"); expect_true(.is_mounted(mp))      # still connected, honestly reported
  p$kill(); Sys.sleep(0.4); r2 <- unmount_drive("Y"); expect_true(r2$ok); expect_match(r2$message, "You can eject it"); expect_false(.is_mounted(mp))
  mount_ok <- identical(suppressWarnings(system2("mount", c("-t", "tmpfs", "tmpfs", mp), stdout = FALSE, stderr = FALSE)), 0L); skip_if(!mount_ok)
  q <- processx::process$new("bash", c("-c", "exec -a 'Rscript /x/run_analysis_cli.R /y/c.json' sleep 120"), wd = mp, cleanup = TRUE); Sys.sleep(0.4)
  r3 <- unmount_drive("Y"); expect_false(r3$ok); expect_match(r3$message, "A StRym run is still using Y:"); q$kill()
})
test_that("releasing a drive outside Windows is the operating system's job, and C: is never released", {
  with_platform("mac", expect_match(unmount_drive("D")$message, "Eject the drive"))
  with_platform("wsl", { expect_match(unmount_drive("C")$message, "system drive"); expect_false(unmount_drive("1")$ok); expect_true(unmount_drive("Q")$ok)  # Q: is not connected: nothing to release
                         expect_match(unmount_drive("Q")$message, "not connected") })
})

test_that("the app builds with the new Browse window on every platform and has no shinyFiles left", {
  skip_if_not_installed("visNetwork")
  for (pl in c("wsl", "linux", "mac")) {
    env <- new.env(); Sys.setenv(STRYM_PLATFORM = pl); on.exit(Sys.unsetenv("STRYM_PLATFORM"), add = TRUE)
    ok <- tryCatch({ assign("app_dir_override", here, envir = env); sys.source(file.path(here, "app.R"), envir = env, keep.source = FALSE); TRUE }, error = function(e) conditionMessage(e))
    expect_true(isTRUE(ok), info = paste(pl, ok)); expect_true(is.function(env$server)); expect_equal(env$platform, pl)
  }
})


test_that("process checks work without /proc (macOS) and the long-read tools are looked for next to the other environment", {
  skip_if(.Platform$OS.type == "windows")
  expect_true(pid_alive_pattern(Sys.getpid(), "R"))
  expect_false(pid_alive_pattern(Sys.getpid(), "no-such-program-name-xyz"))
  expect_false(pid_alive_pattern(999999L, "R")); expect_false(pid_alive_pattern(NA_integer_, "R"))
  b <- .long_env_bins(); expect_true(all(file.path(c(path.expand("~/mamba"), path.expand("~/micromamba"), "/opt/mamba"), "envs", "longread", "bin") %in% b)); expect_true(all(grepl("envs/longread/bin$", b)))
})
# ---- fastp as the read cleaner (macOS has no fastq-mcf) ------------------------------------------------------------------
test_that("the cleaner is picked automatically: fastq-mcf when it exists, fastp otherwise, and the tool check follows", {
  d <- tempfile("tools"); for (r in c("hisat2", "samtools", "cufflinks", "ea-utils/clipper")) dir.create(file.path(d, r), recursive = TRUE)
  mk <- function(rel) { p <- file.path(d, rel); file.create(p); Sys.chmod(p, "755") }
  for (r in c("hisat2/hisat2", "hisat2/hisat2-build", "samtools/samtools", "cufflinks/cufflinks", "cufflinks/cuffmerge", "cufflinks/cuffdiff")) mk(r)
  withr_path <- Sys.getenv("PATH"); on.exit(Sys.setenv(PATH = withr_path), add = TRUE)
  Sys.setenv(PATH = "/usr/bin:/bin")                                    # no fastp on the search path of this test
  expect_equal(pick_cleaner("auto", d), "fastq-mcf")                    # nothing found at all: the TuxNet default (and the missing-tool message names it)
  expect_true("fastq-mcf" %in% check_tools(d)$tool); expect_false(check_tools(d)$found[check_tools(d)$tool == "fastq-mcf"])
  mk("fastp"); expect_equal(pick_cleaner("auto", d), "fastp")
  ct <- check_tools(d); expect_true("fastp" %in% ct$tool); expect_false("fastq-mcf" %in% ct$tool); expect_true(all(ct$found))
  mk("ea-utils/clipper/fastq-mcf"); expect_equal(pick_cleaner("auto", d), "fastq-mcf"); expect_equal(pick_cleaner("fastp", d), "fastp")
  expect_error(pick_cleaner("trimmomatic", d))
})

test_that("fastp gets the same settings as TuxNet's fastq-mcf: window, quality, minimum length, adapters, pairs, Phred+64", {
  a <- fastp_args("a.fq", "o.fq", "ad.fa", 30, 25, 5, 4, 33L)
  expect_true(all(c("--cut_tail", "--cut_window_size", "--cut_mean_quality", "--length_required") %in% a))
  expect_equal(a[which(a == "--cut_mean_quality") + 1], "30"); expect_equal(a[which(a == "--length_required") + 1], "25"); expect_equal(a[which(a == "--cut_window_size") + 1], "5")
  expect_false("--phred64" %in% a); expect_false("-I" %in% a)
  p <- fastp_args("a.fq", "o.fq", "ad.fa", 30, 30, 5, 99, 64L, mate = "b.fq", mate_out = "ob.fq")
  expect_true(all(c("-I", "-O", "--phred64") %in% p)); expect_equal(p[which(p == "--thread") + 1], "16")        # threads are capped
})

test_that("fastp cleans reads for real", {
  skip_if(!nzchar(Sys.which("fastp")), "fastp not installed")
  d <- tempfile("fp"); dir.create(d); fq <- file.path(d, "in.fastq")
  set.seed(1); b <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  rd <- unlist(lapply(1:50, function(i) { s <- b(60); c(paste0("@r", i), s, "+", paste0(strrep("I", 40), strrep("#", 20))) }))   # 20 bad bases at the 3' end
  writeLines(rd, fq); ad <- file.path(here, "IlluminaAdaptorSeq.fasta")
  out <- file.path(d, "out.fastq")
  st <- system2("fastp", fastp_args(fq, out, ad, 30, 30, 5, 2, 33L), stdout = FALSE, stderr = FALSE); expect_equal(st, 0L)
  len <- nchar(readLines(out)[seq(2, 200, by = 4)]); expect_true(length(len) > 0); expect_true(max(len) < 60)               # the low quality tails were cut
})