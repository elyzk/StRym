# syntax=docker/dockerfile:1.7

# This can run on ARM hardware as well thanks to Docker
FROM --platform=linux/amd64 ubuntu:24.04 AS environment

ARG DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/mamba \
    PATH=/root/.local/bin:/opt/mamba/envs/tuxnet/bin:/opt/mamba/envs/longread/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Here we're reusing the predefined dependency installers, we could also just do all of this installation in the Dockerfile
COPY setup_wsl.sh setup_modules.sh /tmp/strym-setup/
RUN bash /tmp/strym-setup/setup_wsl.sh \
    && bash /tmp/strym-setup/setup_modules.sh \
    && micromamba clean --all --yes \
    && rm -rf /tmp/strym-setup /root/.cache


FROM --platform=linux/amd64 ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/mamba \
    HOME=/home/strym \
    STRYM_DATA_ROOT=/data \
    PATH=/opt/mamba/envs/tuxnet/bin:/opt/mamba/envs/longread/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update -qq \
    && apt-get install -y -qq --no-install-recommends ca-certificates curl procps tini \
    && rm -rf /var/lib/apt/lists/* \
    && usermod --login strym --home /home/strym --move-home --shell /bin/bash ubuntu \
    && groupmod --new-name strym ubuntu \
    && install -d -o strym -g strym /data /home/strym/.tuxnet_r

COPY --from=environment /opt/mamba /opt/mamba
COPY . /opt/strym

WORKDIR /opt/strym
USER strym

EXPOSE 3838
STOPSIGNAL SIGTERM

HEALTHCHECK --interval=15s --timeout=5s --start-period=45s --retries=4 \
    CMD /usr/bin/curl --fail --silent --show-error http://127.0.0.1:3838/ >/dev/null || exit 1

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["Rscript", "-e", "app_dir_override <- getwd(); shiny::runApp('app.R', port = 3838, host = '0.0.0.0', launch.browser = FALSE)"]
