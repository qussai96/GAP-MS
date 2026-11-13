FROM mambaorg/micromamba:1.5.8
USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
RUN micromamba install -y -n base -c conda-forge -c bioconda --override-channels --strict-channel-priority python=3.10 pip perl stringtie transdecoder=5.5.0 samtools gffread gffcompare && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH LANG=C.UTF-8 LC_ALL=C.UTF-8
WORKDIR /opt/gapms
COPY pyproject.toml setup.py MANIFEST.in README.md LICENSE THIRD_PARTY.md ./
COPY gapms_v1 ./gapms_v1
RUN python -m pip install --no-cache-dir . && python -c "from importlib.metadata import version; assert version('gapms') == '1.3.2'; assert version('peplocus') == '0.5.6'" && ps -eo pid >/dev/null
WORKDIR /data
ENTRYPOINT ["/opt/conda/bin/gapms"]
LABEL org.opencontainers.image.version="1.3.2" org.opencontainers.image.source="https://github.com/qussai96/GAP-MS"
