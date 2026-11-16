# Simple script improvements we could add later

These are small quality-of-life changes from the StRym v11 scripts. They do not require changing how StRym's packages are divided between environments.

## Use one helper to find StRym's environments

Today, several scripts repeat the same instructions for finding `tuxnet`, `tuxnet-tools`, and `longread`.

A small shared helper could look in the standard StRym locations and then put the right `bin` folders on `PATH`. It should also allow an advanced user to provide a custom environment location.

This would make the launch, test, and command-line scripts shorter and keep them from disagreeing about where StRym is installed.

The helper must support the locations StRym already uses:

- `~/mamba` on Linux and WSL;
- `~/micromamba` on macOS;
- `/opt/mamba` in the Docker image;
- the separate `tuxnet-tools` and `longread` environments.

It should check that R and the required tools are actually present. Finding a folder alone is not enough.

## Give the macOS and Linux installer a clearer name

The current `setup_wsl.sh` also handles macOS, but its name makes that hard to discover.

We could add a small `setup_unix.sh` entry point for macOS and native Linux. It should call the installation logic we already trust rather than introduce a second package list.

`setup_wsl.sh` can remain for compatibility and for Windows users following existing instructions.

## Add a double-click macOS launcher

A short `start_gui.command` file would let a macOS user start StRym from Finder. It would only change into the StRym folder and call `start_gui.sh`.

The first launch may still require right-clicking the file and choosing **Open** because of macOS Gatekeeper.

## Let the launch script open the browser

On native macOS and Linux, `start_gui.sh` could open the StRym page after the server starts. WSL should continue letting `start_gui.bat` open the Windows browser.

The terminal window should plainly say that it must stay open while StRym is running.

## Allow a different port

The launch and stop scripts could read a `STRYM_PORT` setting and use port 3838 when it is not provided.

This would help when another program already uses port 3838 or when two separate StRym instances are needed.

## Make the stop script portable

The stop script could try the tools available on each system in this order:

1. `lsof`, which is normally available on macOS;
2. `ss`, which is common on Linux;
3. `fuser`, as a final Linux fallback.

It should stop only the process listening on StRym's configured port.

## Make long-read tools optional during setup

A setup option such as `--longread` could install minimap2 and featureCounts only for users who need Nanopore or PacBio processing.

Counts-table and short-read users would get a faster, smaller initial installation. The application should continue to give a clear message when a user selects long reads without those tools installed.

## Offer fastp as an optional cleaner

`fastp` can be offered when `fastq-mcf` is unavailable, especially on native macOS. Automatic selection should continue to prefer `fastq-mcf` so Docker and WSL preserve the original TuxNet behavior.

The interface must state that the two cleaners do not produce identical results. A resumed run must not reuse alignments produced by a different cleaner.

## Keep the current package strategy

These script improvements should not copy v11's combined macOS environment or its ARM-Linux package claims.

StRym should keep:

- separate macOS analysis and legacy-tool environments;
- the pinned SAMtools/HTSlib toolchain;
- explicit conda-forge and bioconda channels;
- the existing download retries and archive handling;
- an early, clear rejection of unsupported ARM Linux installations, with amd64 Docker as the supported alternative;
- post-install checks that fail when required software is missing.
