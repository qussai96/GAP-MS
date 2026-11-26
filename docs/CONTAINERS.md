# Apptainer

The published image is a SIF stored on Docker Hub using ORAS:

```bash
apptainer pull gapms.sif oras://docker.io/qussaiab96/gapms:latest
apptainer run --bind /path/to/data:/data gapms.sif -a /data/assembly.fa -p /data/peptides.tsv -g /data/genes.gtf -o /data/results
```

Use the `oras://` address to pull this image. The `latest` tag currently contains
GAPMS 1.4.0. Check the downloaded version with `apptainer run gapms.sif --version`.

The image includes StringTie, TransDecoder 5.5.0 and procps (`ps`) for Nextflow task
metrics. Runs use CPUs; no GPU option is needed.

## Build from source

From the release checkout:

```bash
apptainer build --fakeroot gapms-1.4.0.sif gapms.def
```

Maintainers can publish the SIF with:

```bash
apptainer push gapms-1.4.0.sif oras://docker.io/qussaiab96/gapms:latest
```

A Dockerfile is also provided for users who want to build a Docker/OCI image.
The published ORAS artifact is intended for Apptainer.
