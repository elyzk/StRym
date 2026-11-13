# Moving StRym to another computer

This zip holds **the program** (R code, the browser interface, the analysis steps, the tests, the baseline scripts) and the bundled
*Arabidopsis* TAIR10 reference (`reference/`, ~190 MB). It does **not** hold:

* your sequencing data (FASTQ files) or any project results (the output folders can be many GB: keep them where they are);
* the installed software (R, HISAT2, Cufflinks, Seurat ... about 4–5 GB): the setup script installs it on the new machine;
* StRym's remembered state (last form, job logs): that lives in `~/.tuxnet_r` *inside WSL* and starts fresh on a new machine.

## Windows 10 / 11 (the setup StRym was built and tested on)
1. Install WSL with Ubuntu (once; needs a restart): in PowerShell, `wsl --install -d Ubuntu`. Open Ubuntu once to finish its first-time setup.
2. Unzip this folder anywhere (a path with spaces is fine). Avoid folders that sync to the cloud while you work.
3. Install the software (once, ~20–40 minutes, ~2 GB download). In PowerShell, inside the unzipped folder:
   ```
   wsl -d Ubuntu -u root -e bash setup_wsl.sh
   ```
   `setup_wsl.sh` installs micromamba and a private environment called `tuxnet`. It pins `samtools 1.16.1` on purpose: without the pin,
   conda would pick `samtools 1.5` to fit the newer packages, which changes the Tuxedo pipeline.
   *To reproduce the exact versions used here* (R 4.2.2, samtools 1.16.1, HISAT2 2.2.3, Cufflinks 2.2.1, Seurat 5.1.0, DESeq2 1.38.0 ...)
   use the recorded environment instead of the script's "latest compatible" choice: inside WSL, with `MAMBA_ROOT_PREFIX=/root/mamba`,
   `micromamba create -n tuxnet --file environment-linux-64-explicit.txt` (exact builds) or `-f environment.yml` (exact versions).
4. Start StRym: double-click `start_gui.bat` (it starts the server in WSL and opens http://localhost:3838). Close its window to stop it.
5. Check it works: `wsl -d Ubuntu -u root -e bash tests/run_tests.sh` (about 10–15 minutes; it prints `PASSED: ... FAILED: 0`).
   Single groups of tests: `bash tests/run_tests.sh baselines`, `... pipeline`, `... resume` and so on (the name after `test-`).

An external hard drive plugged in after WSL started must be mounted once: open the **Drives** tab in StRym and press *Mount drive*.

## Linux (not tested, but should work on Ubuntu)
The scripts expect to run as root with the environment in `/root/mamba`. As root: `bash setup_wsl.sh`, then `bash start_gui.sh` and open
http://localhost:3838. The WSL-specific parts (the Drives tab, `/mnt/c` paths) are simply not used.

## Docker on an ARM host
The complete TUX toolchain is available for `linux/amd64`, not `linux/arm64`. Container architecture is fixed when the
container is created; forcing a conda subdirectory inside an existing ARM container is not sufficient. Create an amd64
Ubuntu container from the host shell:

```bash
docker run --platform linux/amd64 -it --name strym-ubuntu -p 3838:3838 \
  -v "$PWD:/opt/StRym" -w /opt/StRym ubuntu:24.04 bash
```

Then, inside the container:

```bash
bash setup_wsl.sh
bash start_gui.sh
```

Open http://localhost:3838 on the host. The named container retains `/root/mamba`; restart it later with
`docker start -ai strym-ubuntu`. Mount data/output directories with additional `-v host-path:container-path` arguments.

`setup_wsl.sh` installs `libarchive-tools` and uses `bsdtar` for the micromamba bootstrap. This avoids GNU tar's
`Function not implemented` failure in some amd64-on-ARM container runtimes.

## macOS
1. Install Apple's command-line developer tools: `xcode-select --install`.
2. On an Apple Silicon Mac, install Rosetta 2: `softwareupdate --install-rosetta`.
3. From the StRym folder, run `bash setup_wsl.sh`. Despite its historical name, the script detects macOS. It creates
   `~/micromamba/envs/tuxnet` for R and `~/micromamba/envs/tuxnet-tools` for the legacy Tuxedo executables. Intel-only
   Bioconda packages run through Rosetta; `fastq-mcf` is built natively from the pinned ea-utils source.
4. Run `bash start_gui.sh`, then open http://localhost:3838. Stop it with `bash stop_gui.sh` or Ctrl-C.
5. Check the installation with `bash tests/run_tests.sh quick`.

## Practical notes
* **Disk space:** about 5 GB for the software, and roughly 22 GB free for a 4-sample project of ~80 million reads (StRym checks before it starts).
* **Time:** with 1 thread, 4 samples of ~80 million reads take roughly 4½ hours up to the end of Cufflinks and Cuffmerge, and Cuffdiff then needs
  more hours (the last real run had not finished Cuffdiff after 2 h 05 min). More threads help the alignment most.
* **Long runs and computer restarts:** a restart (Windows Update, a scheduled reboot) stops WSL and everything in it, including a running analysis.
  The finished parts are kept on disk. Press **Run** again with the same project and StRym *resumes*: it reuses the intact alignments,
  assemblies and merge, and redoes only what is missing (see `DEVELOPMENT_NOTES.md`). Pausing Windows updates for the duration of a long run avoids the interruption.
* **Where to start reading:** `README.md` (what each page does), `analysis/README.md` (the analysis steps), `DEVELOPMENT_NOTES.md` (how it is built).
* **Real timing seen** (2 threads, 7 GB for Ubuntu): 15 single-end samples of ~29 million reads each took 25 h in total (cleaning and alignment ~15 h,
  Cufflinks ~14.5 min per sample, Cuffdiff ~6.5 h). The time estimate StRym shows before a run is known to be too low for large runs.

## Testing and reporting a problem
* **Quick check that an install is healthy:** run `bash tests/run_tests.sh quick` (about 6 minutes: everything except the slow groups).
  From Windows PowerShell, prefix that command with `wsl -d Ubuntu -u root -e`. `bash tests/run_tests.sh edgecases` takes about
  a minute. The whole suite is `bash tests/run_tests.sh` (30–60 minutes, mostly the real Tuxedo tools and single-cell clustering);
  one group at a time works too: `genist`, `pipeline`, `quantify`, `rtpstar`, `seurat`, `resume`.
* **What to send when something goes wrong** (so it can be reproduced): the text of the error on the page; the project's `run_summary.csv`,
  `analysis_run.log` and the `00_design` folder (they record the settings and the sample sheet); and the logs StRym keeps inside Ubuntu in
  `/root/.tuxnet_r/` (`analysis_last_run.log`, `last_run.log`). Never send the FASTQ files or other data; the sample sheet and logs are enough.
* **Single cell data:** one 10x folder that holds several samples is split by the *sample* column of the sample sheet (the design page drafts it from the
  sample labels in the cell barcodes). The Seurat step's log lists every sample it skipped and why, and `cells_per_sample.csv` has one row per sample.
* **Drives (Windows):** the **Drives** tab explains how to refresh, mount, reconnect and release a drive before ejecting it.
* **Reading `.h5` files** needs the `hdf5r` R package, which `setup_wsl.sh` / `setup_modules.sh` now install; on an older install run `setup_modules.sh`.
