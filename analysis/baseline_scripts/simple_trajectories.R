# Simple trajectory plots: Control vs Heat over time, PANA and PI side by side
# Run: Rscript --vanilla ~/simple_trajectories.R
library(ggplot2)

#load folders with data files

folder    <- "/Users/jlreed4/Downloads/genes_read_group_tracking_files"   # tracking + diff files 
genes     <- readLines("/Users/jlreed4/phosphosite_overlap/my_genes.txt")  # one gene ID per line, list of the genes to plot
out_file  <- "/Users/jlreed4/phosphosite_overlap/trajectories.pdf"
drop_reps <- c("Pi_Heated_T1:0", "Pi_Heated_T2:2", "Pana_Control_T3:0",    # PCA outlier replicates
               "Pi_Control_T2:0", "Pana_Control_T4:2", "Pi_Control_T4:2")

genes <- sub("^GLYMA\\.", "Glyma.", trimws(genes[nzchar(trimws(genes))]))  # Cufflinks writes "Glyma."

# read all gene tracking files from cuffdiff; genotype and timepoint come from the file name
files <- list.files(folder, pattern = "read_group_tracking", full.names = TRUE)
d <- do.call(rbind, lapply(files, function(f) {
  x <- read.delim(f)
  x$genotype  <- sub(".*_(Pana|Pi)_T.*", "\\1", basename(f))
  x$file_t    <- as.numeric(sub(".*_T([0-9]+).*", "\\1", basename(f)))
  x$treatment <- ifelse(grepl("heat", x$condition, ignore.case = TRUE), "Heat", "Control")
  x
}))
d$rep_id    <- paste0(d$genotype, ifelse(d$treatment == "Heat", "_Heated_T", "_Control_T"), d$file_t, ":", d$replicate)
d           <- d[!d$rep_id %in% drop_reps, ]
d$timepoint <- paste0("T", d$file_t - 1)              # file T1 = 0 min = T0

# Glyma IDs -> Cufflinks XLOC IDs (via the 'gene' column of one gene_exp.diff file)
diff <- read.delim(list.files(folder, pattern = "gene_exp.diff", full.names = TRUE)[1])
key  <- sapply(genes, function(g) c(diff$test_id[grepl(g, diff$gene, fixed = TRUE)], g)[1])
d    <- d[d$tracking_id %in% key, ]
d$gene <- names(key)[match(d$tracking_id, key)]

# plot: line = mean, points = replicates
p <- ggplot(d, aes(timepoint, FPKM, colour = treatment, group = treatment)) +
  geom_point(alpha = 0.4) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  facet_grid(gene ~ genotype, scales = "free_y") +
  scale_colour_manual(values = c(Control = "grey40", Heat = "firebrick")) +
  expand_limits(y = 0) + theme_bw() +
  theme(strip.text.y = element_text(angle = 0, size = 7))
ggsave(out_file, p, width = 7, height = 1.5 * length(unique(d$gene)) + 1, limitsize = FALSE)
cat("Plotted", length(unique(d$gene)), "of", length(genes), "genes ->", out_file, "\n")
