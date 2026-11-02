# Aho–Corasick protein mapper evaluation

GAP-MS 0.3.0 replaces the former ProteoMapper step with the built-in
`gapms.protein_mapper.ProteinMatcher`. Both prediction and BAM branches use it by
default. The comparisons below were completed on the implementation before the
switch, using GAP-MS 0.2.2 outputs and fresh ProteoMapper runs as references.

## Design

Normalize and deduplicate peptides, canonicalize I/L, and build one C-backed pyahocorasick automaton. Scan each protein sequence independently, retaining all overlapping occurrences and one-based positions. Original query labels are preserved. Unknown residues and stops act as barriers. Queries without hits receive an UNMAPPED row. Protein identifiers use the first FASTA header token; duplicates fail clearly. Sorted output uses GAPMS's existing three-column mapping contract.

One matcher instance is shared between prediction and BAM-derived protein FASTAs. The workflow test asserts one construction and two searches. Different species have different peptide lists and require different automatons. GAP-MS reuses the protein-search automaton; PepLocus genomic mapping still builds its own automaton.

## Experimental protocol

Nine plant species, each with Helixer/FragPipe peptides, were searched against prediction and BAM-derived protein FASTAs: 18 complete comparisons. Compare every (original peptide, protein identifier, one-based position), plus unmapped query sets, against both existing GAPMS ProteoMapper output and fresh ProteoMapper runs. Fresh runs build indexes on private FASTA copies and use eight Perl search threads. The Python mapper uses one process. Large temporary mappings/indexes are removed after comparison; species inputs and outputs remain untouched.

Timings include Python scanning and output separately from setup. Peptide loading plus automaton construction was measured in a separate mapper-only run and is paid once for both branches. Fresh Perl timings include index construction, reading and search/output; copying the input FASTA is excluded. Each configuration was measured once; this is not a replicated scaling study or a benchmark of persistent, already-built ProteoMapper indexes.

## Results

**All 18 comparisons have identical complete match sets and unmapped biological query sets**, against both existing and fresh ProteoMapper outputs.

| Species | One-time Python input/build (s) | Python prediction + BAM scan/output (s) | Perl prediction + BAM index/search (s) | Python peak RSS (MiB) | Perl peak index/search RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| Daucus_carota | 2.74 | 10.81 | 117.27 | 273.7 | 470.1 |
| Spinacia_oleracea | 1.96 | 7.64 | 83.71 | 215.7 | 414.6 |
| Citrus_sinensis | 1.87 | 7.41 | 86.09 | 211.8 | 415.2 |
| Cucumis_melo | 2.16 | 9.30 | 89.36 | 244.7 | 422.9 |
| Solanum_lycopersicum | 2.62 | 8.55 | 102.73 | 243.6 | 431.0 |
| Malus_domestica | 2.09 | 11.15 | 121.74 | 243.0 | 498.3 |
| Zea_mays | 2.76 | 10.87 | 132.32 | 275.2 | 474.0 |
| Hordeum_vulgare | 2.12 | 8.12 | 117.12 | 232.5 | 478.2 |
| Humulus_lupulus | 1.74 | 5.70 | 73.63 | 196.8 | 406.8 |

The combined Python setup-plus-two-search cost is approximately 7.8–11.4 times lower than fresh two-branch Perl cost in these measurements. Setup and scan times come from the separately described measurements, so the ratio is descriptive rather than a replicated end-to-end pipeline speedup. Peak RSS for Python covers the mapper-only process over both searches; Perl reports the larger maximum of its separate indexing/search processes. Reference-comparison set bookkeeping is excluded.

## Downstream equivalence

Carrot was rerun for both branches through feature extraction, classification, fast model metadata attachment and reference-report generation, holding genomic placements, annotation and TMAP fixed. Nineteen files were checked per branch. Feature scores, all protein lists, supported FASTA/GTF, plot images, report TSVs and event-specific FASTA/GTF exports are byte-identical. Prediction metadata BED is byte-identical; BAM metadata has identical records and multiplicities with a different ordering among records tied at the same genomic coordinate.

The only unmapped-list content difference is removal of the literal table header `Peptide Sequence`. ProteoMapper treated it as an invalid unmapped query; the Python input parser correctly consumes it as a header. This does not change biological peptide evidence.

## Tests and scope

The regression suite includes randomized comparisons with an independent
substring-search oracle, I/L aliases, repeated and overlapping occurrences,
stop/unknown-residue boundaries, input normalization, invalid FASTA identifiers,
and one automaton reused by the two workflow branches. Packaging checks verify
that the wheel contains the new mapper and excludes both former Perl scripts.

The evidence supports replacing GAPMS's current exact plain-FASTA mapping step. It does not establish replacement of all ProteoMapper capabilities: PEFF variant annotations and fuzzy/mass-based search options are not implemented. Valid peptides shorter than the default five-AA ProteoMapper key are supported by the new mapper and can differ from legacy results on such inputs. Invalid/ambiguous peptide input is rejected by the shared validated parser. These benchmark lists consist of standard peptide sequences.

## Artifacts and reproduction

Benchmark job 193441; mapper-only memory job 193442; carrot downstream job 193440. The initial jobs 193438/193439 were corrected for shared-path/Perl-environment setup and are not used for performance results. Scripts are in the sibling `protein_mapper_development/` directory: `benchmark.py`, `memory.py`, `downstream.py`, their Slurm scripts, and the tests. `summary.json`, per-species `results.json`, and `carrot_downstream/results.json` hold detailed results.
