# analysis/placeholder.R -- helper used by the placeholder step scripts until the real scripts arrive.
# A placeholder writes a clearly labeled note and picture into the step's output folder, so the whole
# automation (plan → run → organized folder → summary) can be exercised end to end.

.svg_escape <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

strym_placeholder <- function(ctx, id, title, will_do, expects) {
  dir.create(ctx$out_dir, recursive = TRUE, showWarnings = FALSE)
  d <- ctx$design
  note <- c(
    sprintf("PLACEHOLDER — %s", title), "",
    "This step has not been implemented yet: StRym ran a stand-in so the automation can be tested.",
    sprintf("To add the real analysis, put your script in  analysis/%s.R  (see analysis/README.md).", id), "",
    "What this step will do:", paste0("  ", will_do), "",
    "What the real script receives (ctx):",
    sprintf("  data type: %s;  conditions: %d;  time points: %d", d$data_type, d$n_conditions, d$n_timepoints),
    sprintf("  data: %s", if (nzchar(ctx$data_path)) ctx$data_path else "(none given)"),
    sprintf("  sample sheet: %s", if (nzchar(ctx$sample_sheet)) ctx$sample_sheet else "(none given)"),
    if (length(ctx$prereq)) sprintf("  results of earlier steps: %s", paste(names(ctx$prereq), collapse = ", ")),
    "", "Expected outputs:", paste0("  ", expects))
  writeLines(note, file.path(ctx$out_dir, "PLACEHOLDER.txt"))
  lines <- c(.svg_escape(title), "placeholder — real analysis to be added",
             .svg_escape(sprintf("%s · %d conditions · %d time points", if (identical(d$data_type, "bulk")) "whole organ/tissue" else "single cell/nucleus", d$n_conditions, d$n_timepoints)))
  svg <- c('<svg xmlns="http://www.w3.org/2000/svg" width="640" height="360" viewBox="0 0 640 360">',
           '<rect width="640" height="360" fill="#F4F6FC"/><rect width="640" height="64" fill="#4156A1"/>',
           sprintf('<text x="32" y="40" font-family="Arial,sans-serif" font-size="22" fill="#ffffff">%s</text>', lines[1]),
           sprintf('<text x="32" y="150" font-family="Arial,sans-serif" font-size="18" fill="#24347B">%s</text>', lines[2]),
           sprintf('<text x="32" y="190" font-family="Arial,sans-serif" font-size="14" fill="#4D4D4D">%s</text>', lines[3]),
           '<g stroke="#6F7D1C" stroke-width="2" fill="none"><path d="M32 300 C120 240 180 330 260 270 S 420 230 520 290"/></g>',
           '</svg>')
  writeLines(svg, file.path(ctx$out_dir, sprintf("%s_placeholder.svg", id)))
  ctx$log(sprintf("%s: placeholder written (replace analysis/%s.R with the baseline script)", title, id))
  list(status = "placeholder", files = c("PLACEHOLDER.txt", sprintf("%s_placeholder.svg", id)),
       message = "Placeholder — baseline script not added yet.")
}
