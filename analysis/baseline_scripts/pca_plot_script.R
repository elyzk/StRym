#replicate-level PCA plot script (code adapted from Claude, run on the ffar soybean data)
#this script should work with either cuffdiff genes.read_group_tracking output files

#Configuration

#load necessary libraries
library(data.table)
library(ggplot2)

# Read all Cuffdiff replicate FPKM files, this is just T1-T6 for both genotypes
files <- list.files("/Users/jlreed4/Downloads/genes_read_group_tracking_files", "genes\\.read_group_tracking", recursive = TRUE, full.names = TRUE)
d <- rbindlist(lapply(files, fread, select = c("tracking_id", "condition", "replicate", "FPKM"))) #stacking rows of information from the necessary columns contained in each genes.read_group_tracking_file
d[, condition := sub("Dujia", "PI", condition)]
d[, sample := paste0(condition, "_R", replicate + 1)] #create a unique name for each replicate (adding 1 so they start at 1 instead of 0)
d[, .N, by = .(tracking_id, sample)][N > 1] #catch any duplicate samples before assembling the pca

#optional to exclude any outliers 
exclude <- c("T1_PI_Heated_R1", "T2_PI_Heated_R3", "T3_Pana_Control_R1")


#Data Preparation

# Gene x sample matrix; drop excluded; log2(FPKM + 1); 500 most variable genes across all samples to inform the PCA plots
m <- as.matrix(dcast(d, tracking_id ~ sample, value.var = "FPKM", fun.aggregate = sum, fill = 0),
               rownames = "tracking_id") #build matrix
m <- m[, !colnames(m) %in% exclude] #drop excluded samples
m <- log2(m[rowMeans(m) >= 1, ] + 1) #keep samples whose average FPKM across samples is at least 1, compress scale, +1 to avoid log(0)
m <- m[order(apply(m, 1, var), decreasing = TRUE)[1:500], ] #find gene variance across a sample, order those gene variance values from most variance to least, take the first 500


# PCA 

# t(m) to transpose the matrix so genes become columns and samples as rows
#prcomp() treats rows as samples and columns as variables (why we need to transpose m first)
pca <- prcomp(t(m)) 

#square to convert the pca sdev into variance, divide by sum of PCs variance and multiply by 100 to get percentage
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# Labels from condition names like "T6_Pana_Heated"
s <- data.table(sample = rownames(pca$x), pca$x[, 1:2]) #each samples coordinates for the different PCs
s[, cond      := sub("_R[0-9]+$", "", sample)] #removing the replicate suffix so it is just timepoint_genotype_condition
s[, genotype  := ifelse(grepl("pana", cond, ignore.case = TRUE), "Pana", "PI")] #labeling the genotype
s[, treatment := ifelse(grepl("heat", cond, ignore.case = TRUE), "Heat", "Control")] #labeling the treatment
s[, time      := factor(paste0("T", sub(".*T([0-9]+).*", "\\1", cond)), levels = paste0("T", 1:6))] #timepoint

#plotting the PCA
# Circle = Pana, triangle = PI; solid = Heat, hollow = Control; color = timepoint
s[, shp := paste(genotype, treatment)]
ggplot(s, aes(PC1, PC2)) +
  geom_point(aes(colour = time, shape = shp), size = 3, stroke = 1.2) +
  scale_shape_manual(values = c("Pana Control" = 1, "Pana Heat" = 16,
                                "PI Control" = 2, "PI Heat" = 17)) +
  scale_colour_viridis_d(option = "turbo") +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       colour = "Timepoint", shape = NULL) +
  theme_bw()
ggsave("PCA.png", width = 7, height = 5, dpi = 300)

#additional outputs for verification
nrow(s)                 # total points on the plot
table(s$cond)           # replicates per condition
