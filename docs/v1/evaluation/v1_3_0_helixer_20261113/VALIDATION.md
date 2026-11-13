# GAPMS 1.3.0: nine-crop release validation

The release uses PepLocus 0.5.6, a minimal default help screen, and a self-contained pip distribution. The legacy GAPMS wheel and bundled mapper are not runtime dependencies. All nine runs used Helixer FragPipe peptide lists, Helixer annotations and VARUS BAMs, with eight CPUs per run. RefSeq annotations were used only for independent evaluation.

## Results

| Species | Peptides mapped | Ambiguous (%) | Supported models | Specific / ambiguous proposals | Runtime (min) |
|---|---:|---:|---:|---:|---:|
| Daucus carota | 234,952 / 234,952 (100.00%) | 20.7 | 21,052 | 8 / 2,357 | 39.9 |
| Spinacia oleracea | 173,456 / 173,456 (100.00%) | 16.5 | 14,599 | 6 / 1,518 | 34.7 |
| Citrus sinensis | 172,049 / 172,049 (100.00%) | 17.7 | 16,883 | 6 / 1,658 | 33.4 |
| Cucumis melo | 191,103 / 191,103 (100.00%) | 15.3 | 23,098 | 1 / 2,290 | 40.5 |
| Solanum lycopersicum | 200,443 / 200,443 (100.00%) | 18.7 | 18,064 | 4 / 1,721 | 37.3 |
| Malus domestica | 192,384 / 192,384 (100.00%) | 44.3 | 25,278 | 2 / 3,807 | 47.4 |
| Zea mays | 226,092 / 226,092 (100.00%) | 30.5 | 25,551 | 5 / 4,114 | 53.5 |
| Hordeum vulgare | 188,161 / 188,161 (100.00%) | 23.3 | 23,512 | 1 / 2,332 | 28.5 |
| Humulus lupulus | 163,753 / 163,753 (100.00%) | 22.7 | 18,206 | 1 / 2,089 | 18.7 |

Independent Biopython checks verified **5,395,772 peptide placements**, **323,484 peptide-supported RNA ORFs**, and **31,180 diagnostic placement/event links**. All **438,005 projected RNA ORFs** were compared with the original TransDecoder protein predictions. The only differences were 8 residues in 2 unsupported ORFs: deterministic IUPAC codons resolve to specific amino acids in GAPMS/Biopython, while TransDecoder writes X; each difference was independently checked against the genomic codon. All nine audits passed. BED coordinates, blocks, strand labels, all/tryptic multiplicities, supported-annotation structures and score-table schemas were reconciled. Source files retained their original sizes and modification times.

Runtimes include the one-time assessment correction performed during release validation. Runs used shared Slurm nodes with different hardware, so these are descriptive timings rather than a controlled speedup comparison.

## Scientific interpretation

Mapping percentages measure sequence compatibility, not peptide-identification precision. These peptide lists were identified against Helixer proteins, so high recovery is expected and is not an independent demonstration of novel-peptide sensitivity. Most review hypotheses have ambiguous peptide support; this is shown explicitly in the report and its stacked proposal bars. Ambiguity is conditional on the searched assembly and splice paths. Unknown peptide quality is not promoted to high-confidence discovery.

Supported models follow the original high/low rule seeds and iterative classification. They preserve the supplied structures; local peptide support does not validate a complete transcript. The independent RefSeq comparison below uses exact CDS blocks, strand and leading phase, with an unchanged reference denominator. Improved precision after selection can accompany reduced recall. RefSeq is a comparator, not an infallible truth set.

| Species | CDS precision: Helixer → supported (%) | CDS recall: Helixer → supported (%) |
|---|---:|---:|
| Daucus carota | 47.3 → 62.9 | 48.5 → 32.7 |
| Spinacia oleracea | 42.7 → 59.0 | 35.7 → 24.6 |
| Citrus sinensis | 51.8 → 62.4 | 42.1 → 31.2 |
| Cucumis melo | 39.1 → 40.0 | 44.4 → 32.0 |
| Solanum lycopersicum | 37.1 → 64.1 | 43.1 → 31.9 |
| Malus domestica | 50.7 → 63.8 | 45.5 → 32.3 |
| Zea mays | 46.7 → 57.4 | 42.2 → 31.5 |
| Hordeum vulgare | 38.4 → 57.6 | 53.4 → 36.5 |
| Humulus lupulus | 52.4 → 62.7 | 37.5 → 27.3 |

## Fixes verified

- Replaced vendored PepLocus 0.3.1 with the public 0.5.6 API and pinned dependency.
- Preserved annotation junctions and filtered RNA evidence separately, with one RNA junction extraction.
- Kept all exact placements and a separate tryptic BED, with independently computed subset multiplicity.
- Updated one-edit rescue summary ingestion; rescue remains outside exact gene/event support.
- Removed stale build modules and the legacy GAPMS runtime dependency.
- Deduplicated repeated event IDs in public review exports.
- Corrected splice novelty checks to include every supplied model, including other frames and noncoding isoforms; removed 44 false splice-novelty proposals across the nine crops.
- Split specific and ambiguous proposal support in the report, preventing raw candidate totals from implying confirmed discoveries.

Real-data validation exposed a splice-novelty error when a boundary existed in another input model outside the narrow frame-compatible comparison set. The final assessment checks all supplied exon/CDS boundaries on the same chromosome and strand. The assessment stage was rerun from retained verified placements and model occurrences; genomic mapping, RNA reconstruction and classifier fitting were not repeated. Every restored occurrence was checked against its original mapping ID, and non-event model outputs had to reproduce their original hashes. Corrected assessments and report refreshes have separate provenance records in each run manifest. Report refreshes preserved the corrected core scientific files.

## Installation and evidence

```bash
python -m pip install dist/gapms-1.3.0-py3-none-any.whl
gapms --help
```

A fresh pip installation and dependency check passed. The regression suite passed with the automatically resolved runtime. Wheel and source distribution passed Twine checks; rebuilding from the source distribution reproduced the Python modules. The package is prepared locally; it has not been uploaded to PyPI or pushed by this validation run.

See `inputs.json`, `results.tsv`, `runtime_versions.json`, `release_identity.json`, `release_changes.json`, `audit/*/audit.json`, and `jobs.json`. Crop outputs are under `0-raw_files/$crop/outputs/gapms/GAPMS_V1_3_0_helixer_fragpipe`.
