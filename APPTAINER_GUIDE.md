# Apptainer Guide for GAP-MS

## What is an Apptainer Image?
An Apptainer (formerly Singularity) image is a containerized environment that packages your application with all its dependencies, making it portable and reproducible across different systems.

## Steps to Create and Use the Image

### 1. Build the Image

Navigate to the GAP-MS directory and build the image:

**Notes:**
- `gapms.def` is the definition file (recipe) that describes how to build the container
- `gapms.sif` is the output image file (Singularity Image Format)
- Building may take 10-15 minutes depending on your system
- You may need sudo/root privileges or use `--fakeroot` flag if available:

```bash
  apptainer build --fakeroot gapms.sif gapms.def
```



### 2. Use the Image

#### Option A: Run directly (using the runscript)
```bash
apptainer run gapms.sif -g annotations.gtf -f proteins.fasta -p peptides.txt
```

#### Option B: Execute specific commands
```bash
apptainer exec gapms.sif gapms -g annotations.gtf -a assembly.fasta -p peptides.txt
```

#### Option C: Interactive shell
```bash
apptainer shell gapms.sif
# Now you're inside the container
gapms --help
```

### 3. Accessing Your Files

By default, Apptainer mounts your home directory and current working directory. To access files elsewhere:

```bash
apptainer exec --bind /path/to/data:/mnt gapms.sif gapms -g /mnt/annotations.gtf -f /mnt/proteins.fasta -p /mnt/peptides.txt
```

## Understanding the Definition File

The `gapms.def` file has several sections:

- **Bootstrap/From**: Base image (miniconda3 for Python and conda)
- **%files**: Files copied from host to container during build
- **%post**: Commands run during build (installations, setup)
- **%environment**: Environment variables set when container runs
- **%runscript**: Default command when using `apptainer run`
- **%help**: Help text shown with `apptainer run-help gapms.sif`

## Checking the Build

After building, verify the image:

```bash
# Check image info
apptainer inspect gapms.sif

# View help text
apptainer run-help gapms.sif

# Test the command
apptainer exec gapms.sif gapms --help
```

## Sharing the Image

```bash
cd /home/students/q.abbas/Git/GAP-MS && apptainer build --sandbox gapms_sandbox gapms.sif 2>&1 | head -20
```
INFO: Starting build...
INFO: Verifying bootstrap image gapms.sif
INFO: Creating sandbox directory...
INFO: Build complete: gapms_sandbox

```bash
cd /home/students/q.abbas/Git/GAP-MS && apptainer remote login --username qussaiab96 oras://docker.io
```

```bash
cd /home/students/q.abbas/Git/GAP-MS && apptainer push gapms.sif oras://docker.io/qussaiab96/gapms:0.1.3
```