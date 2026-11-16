# Temporary Docker distribution instructions

This file temporarily holds the Docker build, runtime, storage, and distribution documentation.

## Image contents and architecture

The image contains StRym's source, static web assets, `reference/`, `example_data/`, all required R libraries, the legacy short-read tools, and the separate long-read tool environment. Running the image starts the Shiny server on port 3838.

The image is always `linux/amd64` because the complete legacy Tuxedo toolchain is not published for `linux/arm64`. Docker Desktop translates the image on Apple Silicon. The same image runs through Docker on macOS, native Linux, and Windows/WSL.

## Build and run with Compose

From the StRym source directory, create the host data directory and start the service:

```bash
mkdir -p strym-data
docker compose up --build
```

Open http://localhost:3838. Compose binds the server only to the local machine and provides two persistent locations:

| Storage | Container path | Purpose |
|---|---|---|
| Host directory `./strym-data` | `/data` | Inputs, projects, and results directly accessible on the host |
| Named volume `strym-state` | `/home/strym/.tuxnet_r` | Shiny job metadata, saved job configurations, logs, and small files uploaded through the browser |

Lifecycle commands:

```bash
docker compose stop          # stop the container
docker compose start         # restart the existing container
docker compose logs -f       # follow server output
docker compose down          # remove the container; preserve ./strym-data and the state volume
docker compose down -v       # also remove the state volume; ./strym-data is never deleted
```

In StRym, select inputs under `/data` and choose a project output folder under `/data`. A container path such as
`/data/projects/experiment-1` appears on the host as `./strym-data/projects/experiment-1`; users can open, copy, archive,
or share those files normally without `docker cp`. The generated project includes its numbered result directories,
tables, plots, logs, and `index.html`.

The Shiny Results tabs expose selected downloads, such as network tables, but large sequencing outputs should remain in
the mounted data directory rather than being transferred through the browser.

## Path fields inside Docker

StRym's input and output fields work in Docker, but they resolve paths in the container's filesystem namespace. A host
path is usable only after it is mounted, and the value entered in StRym must be the container side of that mount.

| Host path | Docker mount | Path entered in StRym |
|---|---|---|
| `/Users/alice/rna` | `/data/input` | `/data/input` |
| `/home/alice/results` | `/data/output` | `/data/output/project-name` |
| `C:\\Users\\Alice\\RNA` | `/data/input` | `/data/input` |

Do not paste `/Users/...`, `C:\\...`, or an unmounted host `/home/...` path into the containerized app. The built-in
file picker exposes the bind mount as **Data**; select it to browse `/data`. Output written elsewhere in the container,
such as `/home/strym/output`, is part of the container's writable layer and is lost when an `--rm` container is removed.

Keeping stable container paths such as `/data/input` and `/data/output` also makes saved `design.json` files reusable
when the container is recreated on another host.

## What `.tuxnet_r` stores

StRym launches long analyses as child `Rscript` processes. It writes the following application state under `~/.tuxnet_r`:

- `<job>_state.json`: child PID, script, start time, output path, result path, and log path;
- `<job>_config.json`: the submitted design or module configuration;
- `<job>_last_run.log`: analysis progress and final status;
- `uploads/`: small files dropped into the browser when StRym cannot resolve them as mounted paths;
- the corresponding TUX run-state and log files.

Persisting this directory lets a reopened browser reconnect to a job that is still running in the same container, redisplay its progress, restore the most recent Experimental Design form, and retain its logs and browser-uploaded files across container recreation.

The volume does **not** keep a running analysis alive when the container stops. Container shutdown terminates the child processes; persisted PID files are treated as inactive after restart. It also does not contain the project results when those are saved under `/data` or another bind mount.

The state volume is optional for a deliberately stateless container. If it is omitted, StRym creates `~/.tuxnet_r` in the writable container layer; that state disappears when an `--rm` container is removed.

Do not mount one state volume into multiple running StRym containers: job names and state filenames are shared rather than namespaced by container or user.

## Use host input and output directories

Build the image:

```bash
docker build --platform linux/amd64 -t strym:local .
```

Run it with explicit mounts:

```bash
docker run --rm --name strym --platform linux/amd64 \
  -p 127.0.0.1:3838:3838 \
  -v "/absolute/path/to/input:/data/input:ro" \
  -v "/absolute/path/to/output:/data/output" \
  -v strym-state:/home/strym/.tuxnet_r \
  strym:local
```

Select input files under `/data/input` and save projects under `/data/output` in StRym. Input is mounted read-only in this example.

The container runs as the unprivileged user `strym` with UID 1000. Native Linux output directories must be writable by UID 1000. Docker Desktop handles ordinary bind-mount permissions on macOS and Windows/WSL.

Omit the `strym-state` mount if preserving form state, job logs, and browser uploads across container removal is not required.

## Distribute the built image

A registry is the normal distribution mechanism:

```bash
docker tag strym:local registry.example.org/strym:VERSION
docker push registry.example.org/strym:VERSION
```

### Image manifest and host architecture

Every image pushed to a registry already has an OCI/Docker image manifest; Docker creates it during `docker push` or
`docker buildx build --push`. No hand-written manifest file is required.

StRym currently has one real container platform: `linux/amd64`. It runs natively on x86-64 Linux and x86-64
Windows/WSL, and the image has been verified with Docker Desktop's amd64 emulation on Apple Silicon. macOS and Windows
still run a Linux container here, so separate `darwin/*` or `windows/*` images are not required. ARM hosts must select
the amd64 image, as `compose.yaml` already does:

```yaml
platform: linux/amd64
```

Publish the single supported platform directly:

```bash
docker buildx build --platform linux/amd64 \
  --tag registry.example.org/strym:VERSION \
  --push .
```

Users running without the Compose file should select it explicitly:

```bash
docker pull --platform linux/amd64 registry.example.org/strym:VERSION
docker run --platform linux/amd64 registry.example.org/strym:VERSION
```

### ARM Linux

Treat native ARM Linux as a best-effort compatibility target, not a native supported platform. The user needs amd64
QEMU/binfmt emulation only when the Docker host does not already provide it. Docker Desktop bundles cross-architecture
emulation; standalone Docker Engine installations may need host-level registration. Test first:

```bash
docker run --rm --platform linux/amd64 alpine uname -m
```

Success prints `x86_64`. `exec format error` means amd64 emulation is unavailable. An administrator can register it
using the method documented by Docker:

```bash
docker run --privileged --rm tonistiigi/binfmt --install amd64
```

This changes the host's `binfmt_misc` configuration and should not be hidden inside StRym's installation commands.
Even when it works, CPU-heavy StRym analysis is slower under QEMU and remains subject to emulator compatibility.

Docker's guidance: [multi-platform images and QEMU](https://docs.docker.com/build/building/multi-platform/).

### Windows on ARM and WSL

An ARM Windows machine runs an ARM64 WSL 2 kernel, so direct installation inside that WSL distribution has the same
missing-package problem as native ARM Linux. Use Docker Desktop for Windows Arm with its WSL 2 backend instead. Docker
currently labels the Windows Arm release **Early Access**.

Docker Desktop provides cross-architecture emulation inside its managed VM. Do not install QEMU inside the user's
Ubuntu WSL distribution when Docker Desktop owns the Docker daemon; test the Docker Desktop engine with the `alpine`
command above. If Docker Engine was installed directly inside ARM WSL instead, follow the native ARM Linux instructions.

Official requirements: [Docker Desktop on Windows Arm](https://docs.docker.com/desktop/setup/install/windows-install/)
and [Docker Desktop's WSL 2 backend](https://docs.docker.com/desktop/features/wsl/).

A multi-platform manifest containing `linux/amd64` and `linux/arm64` would be useful only after a genuine ARM image
exists. The legacy Bioconda toolchain currently prevents that. An ARM manifest entry must not point to the amd64 image;
that would mislabel the binary architecture rather than make it portable.

Alternatively, export one archive that does not require the source checkout on the recipient's machine:

```bash
docker save strym:local | gzip > strym-linux-amd64.tar.gz
```

Import and run it:

```bash
mkdir -p strym-data
gunzip -c strym-linux-amd64.tar.gz | docker load
docker run --rm --name strym --platform linux/amd64 \
  -p 127.0.0.1:3838:3838 \
  -v "$PWD/strym-data:/data" \
  -v strym-state:/home/strym/.tuxnet_r \
  strym:local
```

## Cytoscape

Cytoscape is optional and is not installed in the container. StRym does not invoke a Cytoscape executable: GENIST and
RTP-STAR write Cytoscape-compatible network tables, while the Shiny Results tabs render the networks directly with the
R `visNetwork` package. Install the Cytoscape desktop application on the host only when its additional layout, styling,
or network-analysis features are needed, then download/import the generated network table.

Bundling the Java desktop GUI into the headless Shiny container would add size and display-server complexity without
supporting StRym's server-side analysis.

## Implementation notes

- `Dockerfile` performs a multi-stage build and copies only the installed Micromamba environments into the runtime stage.
- `tini` runs as PID 1 and reaps the background analysis processes.
- The health check requests `http://127.0.0.1:3838/`.
- The Micromamba bootstrap uses `bsdtar`, avoiding GNU tar's `Function not implemented` failure in some amd64-on-ARM container runtimes.
- `/data` and `~/.tuxnet_r` are writable; the application and bundled reference files under `/opt/strym` are read-only to the runtime user.
