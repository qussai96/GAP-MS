# Implementation record

1. Audit 0.3.0, record commit and run legacy regression tests; isolate a branch
   and worktree.
2. Implement a separate installable namespace and genome-only mode. Validate
   strands, frames, I/L equivalence, multiplicity and coordinates.
3. Add repeated annotation inputs, phased CDS extraction, model projection and
   source-aware model support.
4. Integrate sequential RNA junction extraction, StringTie/TransDecoder and
   peptide-supported ORFs through the same locus evidence store.
5. Add conservative novel-ORF, start, end and splice candidates with diagnostic
   peptides and explicit limitations. Retain unsupported models separately.
6. Test all modes and evidence rules. Keep external reconstruction mocks
   explicit in synthetic workflow tests; validate the real tools on crops.
7. Run nine Casanovo crop datasets through Slurm; inspect completion, counts,
   ambiguity, event evidence and independently verified coordinates. Record
   actual results in VALIDATION.md instead of assuming success from job submission.

Release limits are explicit in SCIENTIFIC_METHODS.md. Whole-transcript phasing,
proteomics FDR estimation, de novo transcript reconstruction without BAM, and
experimental validation of translation initiation remain outside this release.
