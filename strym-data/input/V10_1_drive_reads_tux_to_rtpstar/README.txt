V10_1: reads on a drive, Tux, then RTP-STAR from the Tux results (uses a mounted drive)
=======================================================================================

WHAT THIS TESTS
  - a mounted drive: the whole project (reads, annotation, results) lives on a USB stick / external drive
  - GFF3 annotation instead of GTF (genes.gff3 is the same annotation as genes.gtf)
  - the new step at the end of a Tux run: StRym writes a genes file and a gene expression file for RTP-STAR (tux/rtpstar_inputs/)
  - the TUX results page: download buttons and 'Use them in the RTP-STAR tab'; the RTP-STAR tab: 'Start from your Tux results'

SET UP
  1. Copy this whole folder to the root of a USB stick or external drive. Plug it in.
  2. Open StRym > Drives. Windows: press Refresh; if the drive says 'not mounted', select it and press Mount drive. Mac / Linux: press Refresh.
     The drive must show as 'ready'. Write down the name or letter you see.

RUN 1: Experimental Design page, raw reads
  - Project name: V10_1 drive reads. Output folder: a new folder ON THE DRIVE, e.g. <drive>/V10_1_results.
  - Type of data: whole organ/tissue. What do you have: raw reads (FASTQ). Input data: the 'reads' folder on the drive (8 files, 4 per condition).
  - Annotation: genes.gff3 (NOT the .gtf). Genome: genome.fa. Single-end, short reads. Conditions: 2 (CTL, TRT), time points: 1.
  - Tick DESeq2, Data overview, PCA, Interactive web browser. Run.
  EXPECTED: the pre-run window lists 8 samples in 2 groups, with a time estimate (about 3 minutes); every step ends 'ok';
  01_quantification/tux/rtpstar_inputs/ holds rtpstar_genes_from_tux.csv, rtpstar_expression_from_tux.csv and README_rtpstar_inputs.txt;
  the Quantify line in the results table says 'A genes file and an expression file for RTP-STAR are in tux/rtpstar_inputs/'.
  The expression file has 8 sample columns; the genes file lists the differentially expressed genes (see truth.tsv: genes marked up or down).

RUN 2: the TUX tab on the same reads (optional)
  - Input: the reads folder, output: <drive>/V10_1_tux_tab. After the run the Results tab shows the card 'Inputs for RTP-STAR' with the two files.
  - Press 'Use them in the RTP-STAR tab': the page jumps to RTP-STAR with the genes file and the expression file filled in. Leave the TF file empty.
    EXPECTED: a note that no TF file was given (the run takes longer; it still runs). Run it: a network table named after the genes file
    (rtpstar_genes_from_tux_RTPSTAR_network_for_Cytoscape.txt) and a *_network.png picture appear in the output folder.
  - In the RTP-STAR tab, under 'Start from your Tux results', choose <drive>/V10_1_tux_tab, set 'Smallest change' to 2 and press Build: fewer genes.

UNPLUG TEST (needs the drive)
  - With no run going, press Release drive (Windows) and eject the drive. Plug it in again, press Refresh (and Mount drive on Windows).
  - EXPECTED: the drive is listed again as ready, results are still there. While it is unplugged, typing a path on it in any box gives a plain message.

DEBUGGING HINTS
  - If Quantify fails: read the check list shown before the run, then the run log at the bottom of the page.
  - Is the GFF3 refused? Try genes.gtf; if that works and the GFF3 does not, report it (the file is valid GFF3 with gene / mRNA / exon lines).
