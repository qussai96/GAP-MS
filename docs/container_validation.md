# GAPMS 0.2.2 container validation

Built on 2026-11-01 using Apptainer 1.2.2 and `gapms.def` from source
commit `07a64af`. Both artifacts are local build products, excluded from Git:

- `gapms.sif`: runnable SIF image.
- `gapms_sandbox/`: writable directory image extracted from that SIF.

Commands:

```bash
apptainer build --fakeroot gapms.sif gapms.def
apptainer build --sandbox gapms_sandbox gapms.sif
```

The sandbox conversion's full log was retained before displaying its first 20
lines, avoiding a truncated pipeline that could hide the build exit status.
Build caches and temporary work were placed in `/tmp` where configured.

## Checks

- Definition `%test`: installed GAPMS version 0.2.2, metadata module, pysam and
  pyahocorasick imports, map help, and executable presence for gffread,
  gffcompare, samtools, StringTie and both TransDecoder commands.
- Full regression suite against the SIF's installed package in an isolated
  temporary test directory: **92 passed**. Test dependencies were supplied from
  a temporary directory and are not installed in the production image.
- SIF and sandbox mapping smoke tests covered both strands, a split-codon splice
  peptide, and one-edit rescue. Standard and IGV BED outputs were byte-identical.
- Sandbox package import verified from `/opt/conda/lib/python3.10/site-packages`.

SIF SHA256:

```
941d6f718a8f7d3b19077a82a9930fe6a221a8d7eeeb5d3238384b3efbe49319
```

Local full logs: `container-build-0.2.2.log` and
`container-sandbox-0.2.2.log`. Previous SIF/sandbox backups and temporary build
caches were removed after successful validation to conserve quota.

The GitHub current tree excludes generated carrot test outputs. Historical data
remain in prior commits; Git history was not rewritten. Source regression tests,
tutorial inputs, and compact validation reports remain available.
