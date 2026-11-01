"""Disk-backed locus deduplication, BED12, and auditable result tables."""
from collections import Counter
from contextlib import nullcontext
import csv
import json
import sqlite3

from .sequence import canonical, join_blocks, translate


class HitStore:
    def __init__(self, path, peptides, il_equivalent=True):
        self.peptides = peptides
        self.il_equivalent = il_equivalent
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("CREATE TABLE hits (pid INTEGER, chrom TEXT, strand TEXT, blocks TEXT, matched TEXT, source TEXT, support INTEGER, motifs TEXT, PRIMARY KEY(pid,chrom,strand,blocks)) WITHOUT ROWID")

    def add(self, hits):
        rows = [(pid, chrom, strand, json.dumps(blocks, separators=(",", ":")), matched, source, support, motifs) for pid, chrom, strand, blocks, matched, source, support, motifs in hits]
        self.connection.executemany("""INSERT INTO hits VALUES (?,?,?,?,?,?,?,?) ON CONFLICT(pid,chrom,strand,blocks) DO UPDATE SET
            source=CASE WHEN instr(','||hits.source||',', ','||excluded.source||',')>0 THEN hits.source ELSE hits.source||','||excluded.source END,
            support=max(hits.support,excluded.support),
            motifs=CASE WHEN hits.motifs='' THEN excluded.motifs ELSE hits.motifs END""", rows)
        self.connection.commit()

    def mapped(self):
        return {row[0] for row in self.connection.execute("SELECT DISTINCT pid FROM hits")}

    def finalize(self, directory, reference, table=1, verify=True, limited=None, igv_display=False):
        limited = limited or set()
        counts = Counter(dict(self.connection.execute("SELECT pid,count(*) FROM hits GROUP BY pid")))
        stats = Counter({key: 0 for key in ("mapping_loci", "contiguous_loci", "spliced_loci", "exact_loci", "il_equivalent_loci")})
        with open(directory / "mappings.bed", "w") as bed, open(directory / "mappings.tsv", "w") as tsv, (open(directory / "mappings.igv_display.bed", "w") if igv_display else nullcontext()) as igv:
            if igv is not None:
                igv.write('track name="PepLocus_IGV_display" description="Display only: negative-strand peptide labels reversed; use mappings.bed for analysis" itemRgb="On"\n')
            writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
            writer.writerow(["peptide", "chrom", "start", "end", "strand", "blocks", "matched_sequence", "match_type", "mapping_type", "evidence", "junction_support_reads", "splice_motifs", "loci_for_peptide"])
            query = "SELECT * FROM hits ORDER BY chrom,cast(json_extract(blocks,'$[0][0]') AS INTEGER),strand,pid,blocks"
            for pid, chrom, strand, encoded, matched, source, support, motifs in self.connection.execute(query):
                blocks = json.loads(encoded)
                peptide = self.peptides[pid]
                if verify:
                    observed = translate(join_blocks(reference, chrom, strand, blocks), table)
                    if observed != matched or canonical(observed, self.il_equivalent) != canonical(peptide, self.il_equivalent):
                        raise RuntimeError(f"coordinate verification failed for {peptide} at {chrom}:{blocks}:{strand}")
                start, end = blocks[0][0], blocks[-1][1]
                match_type = "exact" if matched == peptide else "il_equivalent"
                kind = "contiguous" if len(blocks) == 1 else "spliced"
                color = "45,110,180" if kind == "contiguous" else "210,100,35"
                fields = [chrom, start, end, peptide, 1000 if counts[pid] == 1 else 0, strand, start, end, color, len(blocks), ",".join(str(b-a) for a,b in blocks), ",".join(str(a-start) for a,b in blocks)]
                bed.write("\t".join(map(str, fields)) + "\n")
                if igv is not None:
                    display_fields = fields.copy()
                    display_fields[3] = peptide[::-1] if strand == "-" else peptide
                    display_fields[8] = "255,0,0" if kind == "spliced" else "153,153,153"
                    igv.write("\t".join(map(str, display_fields)) + "\n")
                writer.writerow([peptide, chrom, start, end, strand, encoded, matched, match_type, kind, source, support, motifs, counts[pid]])
                stats["mapping_loci"] += 1
                stats[f"{kind}_loci"] += 1
                stats[f"{match_type}_loci"] += 1
        with open(directory / "peptides.tsv", "w") as summary, open(directory / "unmapped.tsv", "w") as unmapped:
            summary.write("peptide\tlength\tloci\tstatus\tnovel_search_limited\n")
            unmapped.write("peptide\tlength\tstatus\tnovel_search_limited\n")
            for pid, peptide in enumerate(self.peptides):
                status = "unmapped" if counts[pid] == 0 else "unique" if counts[pid] == 1 else "multiple"
                summary.write(f"{peptide}\t{len(peptide)}\t{counts[pid]}\t{status}\t{int(pid in limited)}\n")
                if not counts[pid]:
                    unmapped.write(f"{peptide}\t{len(peptide)}\t{status}\t{int(pid in limited)}\n")
        stats.update({"mapped_peptides": len(counts), "unmapped_peptides": len(self.peptides)-len(counts), "uniquely_mapped_peptides": sum(n == 1 for n in counts.values()), "multiply_mapped_peptides": sum(n > 1 for n in counts.values()), "novel_search_limited_peptides": len(limited)})
        return dict(stats)

    def close(self):
        self.connection.close()
