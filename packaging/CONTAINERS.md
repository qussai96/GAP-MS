# Containers

Build from the release checkout, not an older source directory.

```bash
apptainer build --fakeroot gapms-1.3.2.sif gapms.def
apptainer run --bind /path/to/data:/data gapms-1.3.2.sif -a /data/assembly.fa -p /data/peptides.tsv -g /data/genes.gtf -o /data/results
```

Docker Hub accepts Docker/OCI images, not a SIF as a Docker image. To publish:

```bash
docker build -t qussai96/gapms:1.3.2 .
docker run --rm qussai96/gapms:1.3.2 --version
docker login
docker push qussai96/gapms:1.3.2
```

After publication, Apptainer users can also pull the Docker image:

```bash
apptainer pull gapms-1.3.2.sif docker://qussai96/gapms:1.3.2
```

The image includes StringTie, TransDecoder 5.5.0 and procps (`ps`), including
support for Nextflow task metrics. No GPU or `--nv` option is needed.
