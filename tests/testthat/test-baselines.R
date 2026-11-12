# The generalized steps must reproduce the original baseline scripts, and work on the output StRym's Tux step really writes.
run_orig <- function(script, data_dir, out_dir) {
  r <- system2("Rscript", c("--vanilla", shQuote(file.path(here, "analysis", "baseline_scripts", script)), shQuote(data_dir), shQuote(out_dir)), stdout = TRUE, stderr = TRUE)
  st <- attr(r, "status"); if (!is.null(st) && st != 0) stop("original script failed: ", paste(utils::tail(r, 5), collapse = "\n"))
  r
}
same_calls <- function(a, b, label) {                       # compare two result tables gene by gene (the columns both have)
  a <- a[order(a$gene), ]; b <- b[order(b$gene), ]
  expect_equal(a$gene, b$gene, info = paste(label, "genes"))
  for (col in grep("^(call_|log2FC_|q_)", names(a), value = TRUE)) expect_equal(a[[col]], b[[col]], info = paste(label, col))
  for (col in c("n_sig_timepoints", "main_set", "direction")) expect_equal(a[[col]], b[[col]], info = paste(label, col))
}

test_that("modified Shannon entropy reproduces the original mse_simple.R (two genotypes, control drift, transients)", {
  root <- tempfile("syn"); dir.create(root); syn <- synth_tux(root)
  out_o <- file.path(root, "o_mse"); invisible(run_orig("mse_simple.R", syn$orig, out_o))
  p <- strym_default_params()
  data <- strym_load_tux(syn$project); tc <- strym_timecourse(data, p, syn$tps, log = function(...) NULL)
  expect_equal(tc$pairs, data.frame(treated = c("A_Heated", "B_Heated"), control = c("A_Control", "B_Control")))   # genotypes paired by the shared name
  run <- strym_run_mse(tc, p)
  for (g in syn$genotypes) {
    orig <- utils::read.csv(file.path(out_o, sprintf("MSE_%s.csv", g)), stringsAsFactors = FALSE)
    mine <- run$results[[paste0(g, "_Heated")]]
    expect_gt(nrow(orig), 20)                                                                                  # the experiment really has calls to compare
    same_calls(orig, mine, paste("MSE", g))
  }
  om <- utils::read.csv(file.path(out_o, "MSE_main_genes.csv"), stringsAsFactors = FALSE)
  expect_equal(sort(paste(om$genotype, om$gene)), sort(paste(sub("_Heated$", "", sub("_vs_.*", "", run$main$comparison)), run$main$gene)))
})

test_that("Shannon entropy reproduces the original shannon_simple.R", {
  root <- tempfile("syn"); dir.create(root); syn <- synth_tux(root)
  out_o <- file.path(root, "o_sh"); invisible(run_orig("shannon_simple.R", syn$orig, out_o))
  p <- strym_default_params()
  tc <- strym_timecourse(strym_load_tux(syn$project), p, syn$tps, log = function(...) NULL)
  run <- strym_run_shannon(tc, p)
  for (g in syn$genotypes) {
    orig <- utils::read.csv(file.path(out_o, sprintf("shannon_%s.csv", g)), stringsAsFactors = FALSE)
    mine <- run$results[[paste0(g, "_Heated")]]
    expect_gt(nrow(orig), 20); same_calls(orig, mine, paste("Shannon", g))
    expect_equal(sort(orig$H_heat[order(orig$gene)]), sort(mine$H_heat[order(mine$gene)]))
  }
})

test_that("the original behaviour is kept when settings change (baseline off, different cutoffs, dropped replicate)", {
  root <- tempfile("syn"); dir.create(root); syn <- synth_tux(root)
  p <- modifyList(strym_default_params(), list(baseline_first = FALSE, min_log2fc = 1.5, q_cutoff = 0.01))
  tc <- strym_timecourse(strym_load_tux(syn$project), p, syn$tps, log = function(...) NULL)
  r1 <- strym_run_shannon(tc, p)$results[["A_Heated"]]
  expect_true(all(abs(r1$log2FC_T2[r1$call_T2 != 0]) >= 1.5)); expect_true(all(r1$q_T2[r1$call_T2 != 0] < 0.01))
  tc$lfc[[1]][, 1] <- 3; tc$qval[[1]][, 1] <- 0                                                              # pretend the first time point differs strongly
  expect_true(any(strym_run_shannon(tc, p)$results[["A_Heated"]]$call_T0 == 1))                              # not a baseline: it can be called
  expect_true(all(strym_run_shannon(tc, modifyList(p, list(baseline_first = TRUE)))$results[["A_Heated"]]$call_T0 == 0))   # a baseline: never called
  d <- strym_load_tux(syn$project); drop <- d$samples$sample[d$samples$condition == "A_Heated" & d$samples$time_point == "T2"][1]
  d2 <- strym_drop_samples(d, list(drop_reps = drop)); expect_false(drop %in% d2$samples$sample); expect_false(drop %in% d2$rep$sample)
  expect_warning(strym_drop_samples(d, list(drop_reps = "not_a_sample")), "ignored")
})

test_that("the pairs of treated and control conditions are found, or the reason is explained", {
  S <- data.frame(condition = c("Ctrl", "Heat", "Cold"), stringsAsFactors = FALSE)
  pr <- strym_pairs(S, strym_default_params()); expect_equal(pr, data.frame(treated = c("Heat", "Cold"), control = c("Ctrl", "Ctrl")))
  expect_equal(strym_pairs(S, modifyList(strym_default_params(), list(control = "Heat")))$control, c("Heat", "Heat"))
  expect_equal(strym_pairs(S, modifyList(strym_default_params(), list(pairs = "Cold:Heat")))$treated, "Cold")
  expect_error(strym_pairs(data.frame(condition = c("B", "D")), strym_default_params()), "could not tell which condition is the control")
  expect_error(strym_pairs(S, modifyList(strym_default_params(), list(control = "Nope"))), "not one of the conditions")
  expect_error(strym_pairs(S, modifyList(strym_default_params(), list(pairs = "Heat"))), "treated:control")
  amb <- data.frame(condition = c("WT_x", "Mock_y", "Heat", "Cold"))
  expect_error(strym_pairs(amb, modifyList(strym_default_params(), list(control_pattern = "WT|Mock"))), "cannot tell which treated condition goes with which")
})

test_that("settings are read back safely from JSON and bad values fall back to the defaults", {
  p <- strym_normalize_params(list(q_cutoff = "0.01", baseline_first = "false", sc_integration = "weird", min_sig_tp = "abc", drop_reps = "a, b"))
  expect_equal(p$q_cutoff, 0.01); expect_false(p$baseline_first); expect_equal(p$sc_integration, "cca"); expect_equal(p$min_sig_tp, 2); expect_equal(p$drop_reps, "a, b")
  expect_equal(names(strym_normalize_params(NULL)), names(strym_default_params()))
  d <- strym_normalize_design(list(params = list(q_cutoff = 0.2))); expect_equal(d$params$q_cutoff, 0.2); expect_equal(d$params$s_max, 3)
})

test_that("time points are put in the order of the experiment", {
  expect_equal(strym_order_timepoints(c("4h", "0h", "1h")), c("0h", "1h", "4h"))
  expect_equal(strym_order_timepoints(c("late", "early"), c("early", "late")), c("early", "late"))
  expect_equal(strym_order_timepoints(c("b", "a")), c("b", "a"))
  expect_equal(strym_order_timepoints(c("T10", "T2", "T1")), c("T1", "T2", "T10"))
})

test_that("modified Shannon entropy copes with few time points", {
  m <- matrix(c(1, 1.1, 9, 0.9, 2, 2, 2, 2), 2, 4, byrow = TRUE, dimnames = list(c("g1", "flat"), paste0("T", 1:4)))
  o <- strym_mse_outliers(m); expect_equal(rownames(o), "g1"); expect_equal(unname(o["g1", ]), c(0, 0, 1, 0))     # the flat gene is dropped, the spike found
  expect_true(all(strym_mse_outliers(m[, 1:2]) == 0))                                                           # 2 points: nothing can be an outlier
  expect_equal(unname(strym_mse_outliers(m[1, c(1, 3, 4), drop = FALSE])["g1", ]) %in% c(-1, 0, 1), rep(TRUE, 3)) # 3 points: runs without error
  expect_true(all(strym_mse_outliers(matrix(0, 3, 5, dimnames = list(letters[1:3], NULL))) == 0))
})

test_that("the comparisons are checked from the sample sheet before a run", {
  sheet <- data.frame(condition = rep(c("Control", "Heat"), each = 6), time_point = rep(rep(c("0h", "1h", "4h"), each = 2), 2))
  ok <- strym_check_comparisons(sheet, strym_default_params(), c("0h", "1h", "4h"))
  expect_equal(ok$pairs, data.frame(treated = "Heat", control = "Control")); expect_length(ok$problems, 0)
  gap <- sheet[!(sheet$condition == "Heat" & sheet$time_point == "4h"), ]
  expect_match(strym_check_comparisons(gap, strym_default_params())$problems, "Condition “Heat” has no samples at time point 4h")
  one <- strym_check_comparisons(sheet[sheet$time_point == "0h", ], strym_default_params()); expect_match(one$problems, "needs? at least 2|need at least 2")
  nc <- strym_check_comparisons(data.frame(condition = c("B", "D"), time_point = "t1"), strym_default_params())
  expect_null(nc$pairs); expect_true(any(grepl("could not tell which condition is the control", nc$problems)))
  expect_equal(strym_check_comparisons(NULL, strym_default_params()), list(pairs = NULL, problems = character()))
})
