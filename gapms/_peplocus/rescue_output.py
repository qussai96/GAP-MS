"""Separate, disk-deduplicated one-edit results; exact outputs stay unchanged."""
from collections import Counter
from contextlib import nullcontext
import csv
import json
import sqlite3

from .rescue import one_edit
from .sequence import join_blocks, translate


class RescueStore:
    def __init__(self, path, peptides, unmapped_ids, eligible_ids, il_equivalent=True):
        self.peptides, self.unmapped_ids = peptides, tuple(unmapped_ids)
        self.eligible = set(eligible_ids)
        self.il_equivalent = il_equivalent
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("""CREATE TABLE hits (
          pid INTEGER, chrom TEXT, strand TEXT, blocks TEXT, matched TEXT,
          source TEXT, support INTEGER, motifs TEXT, kind TEXT, edit_pos INTEGER,
          query_aa TEXT, reference_aa TEXT, cigar TEXT,
          PRIMARY KEY(pid,chrom,strand,blocks)) WITHOUT ROWID""")

    def add(self, hits):
        rows = []
        for pid,chrom,strand,blocks,matched,source,support,motifs,edit in hits:
            rows.append((pid,chrom,strand,json.dumps(blocks,separators=(",",":")),matched,
                         source,support,motifs,edit["edit_type"],edit["edit_query_position"],
                         edit["query_aa"],edit["reference_aa"],edit["alignment_cigar"]))
        self.connection.executemany("""INSERT INTO hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
          ON CONFLICT(pid,chrom,strand,blocks) DO UPDATE SET
          source=CASE WHEN instr(','||hits.source||',',','||excluded.source||',')>0
            THEN hits.source ELSE hits.source||','||excluded.source END,
          support=max(hits.support,excluded.support)""",rows)
        self.connection.commit()

    def finalize(self, directory, reference, table=1, verify=True, igv_display=False):
        counts = dict(self.connection.execute("SELECT pid,count(*) FROM hits GROUP BY pid"))
        stats = Counter({key:0 for key in ("mapping_loci","contiguous_loci","spliced_loci",
                        "substitution_loci","insertion_loci","deletion_loci")})
        with open(directory/"rescue_mappings.bed","w") as bed, open(directory/"rescue_mappings.tsv","w") as tsv, (open(directory/"rescue_mappings.igv_display.bed","w") if igv_display else nullcontext()) as visual:
            if visual:
                visual.write('track name="PepLocus_one_edit_IGV_display" description="One-edit rescue; negative-strand labels reversed for display only" itemRgb="On"\n')
            writer=csv.writer(tsv,delimiter="\t",lineterminator="\n")
            writer.writerow(["peptide","chrom","start","end","strand","blocks","matched_sequence",
                             "match_type","mapping_type","evidence","junction_support_reads","splice_motifs",
                             "loci_for_peptide","edit_distance","edit_type","edit_query_position","query_aa",
                             "reference_aa","alignment_cigar"])
            query="SELECT * FROM hits ORDER BY chrom,cast(json_extract(blocks,'$[0][0]') AS INTEGER),strand,pid,blocks"
            for pid,chrom,strand,encoded,matched,source,support,motifs,kind,position,qaa,raa,cigar in self.connection.execute(query):
                peptide=self.peptides[pid];blocks=json.loads(encoded)
                if verify:
                    observed=translate(join_blocks(reference,chrom,strand,blocks),table)
                    edit=one_edit(peptide,observed,self.il_equivalent)
                    if observed!=matched or edit is None or edit["alignment_cigar"]!=cigar:
                        raise RuntimeError(f"one-edit coordinate verification failed for {peptide} at {chrom}:{blocks}:{strand}")
                start,end=blocks[0][0],blocks[-1][1]
                category="contiguous" if len(blocks)==1 else "spliced"
                fields=[chrom,start,end,peptide,1000 if counts[pid]==1 else 0,strand,start,end,"190,125,25",len(blocks),
                        ",".join(str(b-a) for a,b in blocks),",".join(str(a-start) for a,b in blocks)]
                bed.write("\t".join(map(str,fields))+"\n")
                if visual:
                    f=fields.copy();f[3]=peptide[::-1] if strand=="-" else peptide
                    f[8]="255,0,0" if category=="spliced" else "153,153,153"
                    visual.write("\t".join(map(str,f))+"\n")
                writer.writerow([peptide,chrom,start,end,strand,encoded,matched,"one_edit",category,source,support,
                                 motifs,counts[pid],1,kind,position,qaa,raa,cigar])
                stats["mapping_loci"]+=1;stats[category+"_loci"]+=1;stats[kind+"_loci"]+=1
        with open(directory/"rescue_peptides.tsv","w") as handle:
            writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
            writer.writerow(["peptide","length","eligible","rescue_loci","status"])
            for pid in self.unmapped_ids:
                n=counts.get(pid,0)
                status="not_eligible_short" if pid not in self.eligible else "unmapped" if n==0 else "unique" if n==1 else "multiple"
                writer.writerow([self.peptides[pid],len(self.peptides[pid]),int(pid in self.eligible),n,status])
        stats.update({"exact_unmapped_peptides":len(self.unmapped_ids),"eligible_peptides":len(self.eligible),
                      "short_peptides_not_searched":len(self.unmapped_ids)-len(self.eligible),
                      "rescued_peptides":len(counts),"unresolved_eligible_peptides":len(self.eligible)-len(counts),
                      "uniquely_rescued_peptides":sum(n==1 for n in counts.values()),
                      "multiply_rescued_peptides":sum(n>1 for n in counts.values())})
        return dict(stats)

    def close(self):
        self.connection.close()
