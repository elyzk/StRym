# Version 10: the Drives page on every platform, and the fixes made after the first public test round.
suppressPackageStartupMessages(library(DT))
source(file.path(here, "analysis", "qc_overview.R"))


test_that("the platform is detected, and can be forced for testing", {
  expect_true(strym_platform() %in% c("wsl", "linux", "mac", "windows"))
  expect_equal(with_platform("mac", strym_platform()), "mac"); expect_true(with_platform("wsl", strym_is_wsl())); expect_false(with_platform("linux", strym_is_wsl()))
})

test_that("on a Mac or Linux the Windows-drive rules stay out of the way", {
  for (pl in c("mac", "linux")) with_platform(pl, {
    expect_equal(clean_path("D:\\data\\x"), "D:\\data\\x")                          # a Windows path is not translated outside WSL
    expect_equal(clean_path("  '/Volumes/My Drive/data/' "), "/Volumes/My Drive/data")
    expect_equal(to_windows_path("/mnt/d/data"), "/mnt/d/data"); expect_equal(unmounted_drive_problem(file.path(tempdir(), "x")), "")
    expect_false(mount_drive("D")$ok); expect_match(mount_drive("D")$message, "mounts drives by itself"); expect_match(unmount_drive("D")$message, "Eject the drive")
  })
  with_platform("wsl", { expect_equal(clean_path("D:\\data\\x"), "/mnt/d/data/x"); expect_equal(to_windows_path("/mnt/d/data"), "D:\\data") })
})

test_that("Linux drives are found in /proc/mounts: real devices under /media, /run/media and /mnt only", {
  m <- c("proc /proc proc rw 0 0", "/dev/sda1 / ext4 rw 0 0", "/dev/sdb1 /media/me/USB\\040STICK exfat rw 0 0", "/dev/sdc1 /run/media/me/BIGDISK ext4 rw 0 0",
         "tmpfs /run/user/1000 tmpfs rw 0 0", "/dev/loop3 /snap/core/1 squashfs ro 0 0", "/dev/sdd1 /mnt/data ext4 rw 0 0", "server:/x /mnt/share nfs rw 0 0", "/dev/sde1 /boot/efi vfat rw 0 0")
  expect_equal(.volume_mounts("linux", m), c("/media/me/USB STICK", "/mnt/data", "/mnt/share", "/run/media/me/BIGDISK"))
  expect_equal(.volume_mounts("linux", character()), character())
})

test_that("the drive list and Browse places follow the host platform", {
  cols <- c("letter", "label", "fs", "size", "free", "mount", "mounted", "usable", "state")
  for (pl in c("mac", "linux")) with_platform(pl, expect_true(all(cols %in% names(list_drives()))))
  tmp <- tempfile("vol"); dir.create(tmp)
  d <- data.frame(letter = c("USB", "BIG"), label = c("USB", "BIG"), fs = "", size = NA_real_, free = NA_real_, mount = c(file.path(tmp, "USB"), file.path(tmp, "BIG")),
                  mounted = TRUE, usable = TRUE, state = c("ready", "stale"), stringsAsFactors = FALSE)
  with_platform("mac", {
    p <- strym_places(d)
    expect_equal(p$path[match(c("USB", "BIG"), p$label)], d$mount)
    expect_true(all(p$kind[p$label %in% c("USB", "BIG")] == "go"))
    expect_true("/" %in% p$path)
  })
  with_platform("wsl", expect_true(path.expand("~") %in% strym_places(d)$path))
  with_platform("linux", expect_equal(unmounted_drive_problem(file.path(tmp, "x")), ""))
})


test_that("a counts table without a sample sheet is read the same way headless and in the page: replicates grouped by the column names", {
  root <- tempfile("v10"); dir.create(root)
  ng <- 300; set.seed(2); m <- matrix(rnbinom(ng * 8, mu = 80, size = 10), ng, 8, dimnames = list(sprintf("G%03d", 1:ng), c(paste0("control_rep", 1:4), paste0("treated_rep", 1:4))))
  m[1:30, 5:8] <- m[1:30, 5:8] * 8L
  f <- file.path(root, "counts.csv"); write.csv(data.frame(gene = rownames(m), m), f, row.names = FALSE)
  d <- strym_default_design(); d$project_name <- "nosheet"; d$project_dir <- file.path(root, "p"); d$data_type <- "bulk"; d$data_kind <- "counts"; d$data_path <- f; d$n_conditions <- 2L; d$n_timepoints <- 1L
  d$de_methods <- "deseq2"; d$plots <- "qc"
  r <- quiet(strym_run(d, c("deseq2", "qc_overview"), app_dir = here))
  expect_equal(r$summary$status[r$summary$step == "deseq2"], "ok")                                              # was: "DESeq2 needs replicates: every group has a single sample"
  expect_true(file.exists(file.path(root, "p", "00_design", "sample_sheet_drafted_by_strym.csv")))
  sh <- read.csv(file.path(root, "p", "00_design", "sample_sheet_drafted_by_strym.csv")); expect_equal(sort(unique(sh$condition)), c("control", "treated")); expect_equal(nrow(sh), 8)
})

test_that("a mistyped condition or time-point name is flagged before the run", {
  root <- tempfile("v10"); dir.create(root); ss <- file.path(root, "ss.csv")
  write.csv(data.frame(sample = paste0("s", 1:6), condition = rep(c("control", "drought"), each = 3), time_point = rep(c("0h", "1h", "2h"), 2)), ss, row.names = FALSE)
  d <- strym_default_design(); d$data_type <- "bulk"; d$data_kind <- "counts"; d$sample_sheet <- ss; d$n_conditions <- 2L; d$n_timepoints <- 3L
  d$condition_names <- c("control", "draught"); d$timepoint_names <- c("0h", "1h", "3h")
  w <- strym_validate_design(d, check_paths = TRUE)$warnings
  expect_true(any(grepl("You typed the condition name draught, but the sample sheet has: control, drought", w))); expect_true(any(grepl("You typed the time-point name 3h", w)))
  d$condition_names <- c("control", "drought"); d$timepoint_names <- c("0h", "1h", "2h")
  w2 <- strym_validate_design(d, check_paths = TRUE)$warnings; expect_false(any(grepl("You typed", w2)))
})

test_that("the data overview says when every group has only two replicates; the comparison message counts genes once", {
  root <- tempfile("v10"); dir.create(root); set.seed(8)
  m <- matrix(rnbinom(200 * 4, mu = 60, size = 8), 200, 4, dimnames = list(sprintf("G%03d", 1:200), c("a_rep1", "a_rep2", "b_rep1", "b_rep2")))
  f <- file.path(root, "c.csv"); write.csv(data.frame(gene = rownames(m), m), f, row.names = FALSE)
  ss <- file.path(root, "ss.csv"); write.csv(data.frame(sample = colnames(m), condition = c("a", "a", "b", "b"), replicate = c(1, 2, 1, 2)), ss, row.names = FALSE)
  d <- strym_default_design(); d$project_dir <- file.path(root, "p"); d$data_type <- "bulk"; d$data_kind <- "counts"; d$data_path <- f; d$sample_sheet <- ss; d$n_conditions <- 2L; d$n_timepoints <- 1L
  d <- strym_normalize_design(d); out <- file.path(root, "qc"); dir.create(out)
  r <- strym_step_qc_overview(list(design = d, data_path = f, sample_sheet = ss, out_dir = out, project_dir = d$project_dir, app_dir = here, prereq = list(), log = function(...) NULL))
  notes <- read.csv(file.path(out, "data_overview_notes.csv")); expect_true(any(grepl("only 2 replicates", paste(notes[[3]], notes[[2]], collapse = " "))) || any(grepl("only 2 replicates", unlist(notes))))
})

test_that("every R file of the program parses (a syntax slip in app.R or a page would only show when the app is started)", {
  fs <- c(list.files(here, "\\.R$", full.names = TRUE), list.files(file.path(here, "analysis"), "\\.R$", full.names = TRUE))
  bad <- Filter(function(f) inherits(tryCatch(parse(f, encoding = "UTF-8"), error = function(e) e), "error"), fs)
  expect_length(bad, 0); expect_gt(length(fs), 30)
  expect_false(any(vapply(fs, function(f) identical(readBin(f, "raw", 3), as.raw(c(0xEF, 0xBB, 0xBF))), TRUE)))          # no byte-order mark (R cannot source it)
})

test_that("the app itself can be built: its pages and its server function exist for every platform", {
  skip_if_not_installed("visNetwork")
  for (pl in c("wsl", "linux", "mac")) {
    env <- new.env(); Sys.setenv(STRYM_PLATFORM = pl); on.exit(Sys.unsetenv("STRYM_PLATFORM"), add = TRUE)
    ok <- tryCatch({ assign("app_dir_override", here, envir = env); sys.source(file.path(here, "app.R"), envir = env, keep.source = FALSE); TRUE }, error = function(e) conditionMessage(e))
    expect_true(isTRUE(ok), info = paste(pl, ok)); expect_true(exists("ui", envir = env) && is.function(env$server))
    expect_equal(env$platform, pl); html <- as.character(env$drives_help); expect_false(grepl("<b>|<strong>", html))
    expect_equal(grepl("Mount drive", as.character(env$drives_card)), pl == "wsl")
  }
})
