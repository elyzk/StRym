# Drives: what "mounted" really means, reconnecting a drive that was unplugged and plugged back in, releasing a drive before an eject,
# and messages a biologist can act on. Windows and the drives are replaced by stand-ins, so this runs anywhere.

with_mocks <- function(mocks, code) {          # replace functions in the global environment for the length of `code`, then put everything back
  had <- vapply(names(mocks), function(n) exists(n, envir = globalenv(), inherits = FALSE), logical(1))
  old <- lapply(names(mocks), function(n) if (exists(n, envir = globalenv(), inherits = FALSE)) get(n, envir = globalenv()) else NULL)
  for (n in names(mocks)) assign(n, mocks[[n]], envir = globalenv())
  on.exit(for (i in seq_along(old)) if (had[i]) assign(names(mocks)[i], old[[i]], envir = globalenv()) else rm(list = names(mocks)[i], envir = globalenv()))
  force(code)
}
fake_windows <- function() data.frame(letter = c("C", "D", "E"), label = c("Windows", "LaCie", ""), fs = "NTFS", size = c(5e11, 2e12, 1e12), free = c(6e10, 1.6e12, 9e11), stringsAsFactors = FALSE)

test_that("a mounted drive only counts when it answers; the answer can be remembered briefly", {
  expect_true(.mount_usable(tempdir())); expect_false(.mount_usable(file.path(tempdir(), "no_such_drive")))
  d <- file.path(tempdir(), "later"); expect_false(.mount_usable(d, ttl = 60)); dir.create(d)
  expect_false(.mount_usable(d, ttl = 60))                                        # remembered
  expect_true(.mount_usable(d))                                                  # asked again
  .forget_usable(); expect_true(.mount_usable(d, ttl = 60))
})

test_that("every drive gets a state: ready, not mounted, stale (unplugged and replugged) or removed", {
  with_mocks(list(.windows_volumes = function() fake_windows(), dir.create = function(...) TRUE,
                  .is_mounted = function(mp) mp %in% c("/mnt/c", "/mnt/d", "/mnt/f"),
                  .wsl_drive_mounts = function() c("/mnt/c", "/mnt/d", "/mnt/f"),
                  .mount_usable = function(mp, ...) mp == "/mnt/c"), {
    d <- list_drives(); st <- setNames(d$state, d$letter)
    expect_equal(st[["C"]], "ready"); expect_equal(st[["D"]], "stale"); expect_equal(st[["E"]], "not mounted"); expect_equal(st[["F"]], "removed")
    expect_equal(d$usable[d$letter == "C"], TRUE); expect_equal(d$usable[d$letter == "D"], FALSE); expect_equal(d$label[d$letter == "F"], "")
    r <- drive_roots(d)                                                           # a drive Windows no longer has is not offered in Browse…
    expect_true(all(c("Home", "Root") %in% names(r))); expect_true(all(c("/mnt/c", "/mnt/d", "/mnt/e") %in% r)); expect_false("/mnt/f" %in% r)
    expect_true("D: LaCie" %in% names(r))
  })
  with_mocks(list(.windows_volumes = function() NULL), {                         # Windows cannot be asked: no drive is called "removed"
    expect_false("removed" %in% list_drives()$state)
  })
})

test_that("mount_drive reconnects a dead connection, explains failures, and leaves a working drive alone", {
  expect_false(mount_drive("1")$ok); expect_match(mount_drive("")$message, "Pick a drive letter")
  st <- new.env(); st$mounted <- TRUE; st$dead <- TRUE; st$calls <- list()
  fake <- list(.is_mounted = function(mp) st$mounted, .mount_usable = function(mp, ...) !st$dead, dir.create = function(...) TRUE,
               system2 = function(cmd, args = character(), ...) {
                 st$calls[[length(st$calls) + 1]] <- c(cmd, args)
                 if (cmd == "umount") st$mounted <- FALSE
                 if (cmd == "mount") { st$mounted <- TRUE; st$dead <- FALSE }
                 character() })
  with_mocks(fake, {
    r <- mount_drive("d"); expect_true(r$ok); expect_match(r$message, "^Mounted D: at /mnt/d"); expect_match(r$message, "D:\\\\data")
    expect_equal(vapply(st$calls, `[`, "", 1), c("umount", "mount")); expect_equal(st$calls[[1]][2], "-l")           # the dead connection is dropped first
    expect_equal(st$calls[[2]][2:4], c("-t", "drvfs", "D:"))
    st$calls <- list(); again <- mount_drive("D"); expect_true(again$ok); expect_match(again$message, "already connected and working"); expect_length(st$calls, 0)
  })
  st$mounted <- FALSE; st$dead <- TRUE; st$calls <- list()
  with_mocks(modifyList(fake, list(system2 = function(cmd, args = character(), ...) "mount: /mnt/z: special device Z: does not exist."), keep.null = TRUE), {
    r <- mount_drive("Z"); expect_false(r$ok); expect_match(r$message, "Windows has no drive Z: right now"); expect_match(r$message, "plugged in"); expect_match(r$message, "Refresh")
  })
})

test_that("unmount_drive lets go of a drive for ejecting, refuses C:, and detaches a busy one", {
  expect_false(unmount_drive("C")$ok); expect_match(unmount_drive("C")$message, "system drive"); expect_false(unmount_drive("1")$ok)
  st <- new.env(); st$mounted <- FALSE; st$busy <- FALSE; st$calls <- list()
  fake <- list(.is_mounted = function(mp) st$mounted,
               system2 = function(cmd, args = character(), ...) {
                 st$calls[[length(st$calls) + 1]] <- c(cmd, args)
                 if (cmd == "umount" && (!st$busy || "-l" %in% args)) st$mounted <- FALSE
                 if (st$busy && !"-l" %in% args) "umount: /mnt/d: target is busy." else character() })
  with_mocks(fake, {
    r <- unmount_drive("D"); expect_true(r$ok); expect_match(r$message, "not connected"); expect_length(st$calls, 0)           # nothing to release
    st$mounted <- TRUE; r <- unmount_drive("d"); expect_true(r$ok); expect_match(r$message, "Released D:.*eject it in Windows"); expect_length(st$calls, 1)
    st$mounted <- TRUE; st$busy <- TRUE; st$calls <- list(); r <- unmount_drive("D"); expect_true(r$ok); expect_length(st$calls, 2); expect_equal(st$calls[[2]][2], "-l")
  })
})

test_that("a path on a drive that is not connected, or connected but dead, gets a plain message with the fix", {
  with_mocks(list(.is_mounted = function(mp) mp == "/mnt/d", .mount_usable = function(mp, ...) FALSE), {
    expect_match(unmounted_drive_problem("/mnt/e/data"), "Drive E: is not connected to StRym yet.*Refresh.*Mount drive")
    expect_match(unmounted_drive_problem("/mnt/d/data/x.csv"), "Drive D: is connected but not responding.*unplugged or ejected.*Mount drive")
    expect_equal(unmounted_drive_problem("/root/data"), ""); expect_equal(unmounted_drive_problem("/mnt/wsl/x"), "")
  })
  with_mocks(list(.is_mounted = function(mp) TRUE, .mount_usable = function(mp, ...) TRUE), expect_equal(unmounted_drive_problem("/mnt/d/data"), ""))
})

test_that("on this computer the real drive list is sane (skipped outside WSL)", {
  skip_if(!.is_mounted("/mnt/c"), "no /mnt/c here")
  d <- list_drives(); expect_true(all(c("letter", "mounted", "usable", "state") %in% names(d)))
  if (nrow(d) && "C" %in% d$letter) expect_equal(d$state[d$letter == "C"], "ready")
  expect_true(mount_drive("C")$ok)
  with_mocks(list(dir.create = function(...) TRUE), { q <- mount_drive("Q"); expect_false(q$ok); expect_match(q$message, "Windows has no drive Q|Could not mount Q") })   # no drive Q: the message says so
  expect_equal(unmounted_drive_problem("/mnt/c/Users"), ""); expect_true(all(c("Home", "Root") %in% names(drive_roots())))
})
