"""Fast attachment of protein coordinates to PepLocus-compatible genomic loci.

The assembly mapper remains authoritative. This module preserves the historical
metadata BED contract used by reference reports without translating/projecting
an entire protein afresh for every peptide.
"""
from bisect import bisect_right
from pathlib import Path
import time

import numpy as np
from Bio import SeqIO

from .gtf_utils import gtf_to_df_with_genes
from .peptide_utils import mapping_file_to_df, get_gene_protein_specific_peps, resolve_peptide_starts
from .protein_projection import (_coding_nt_positions, classify_peptide_type,
                                 export_to_bed, get_protein_to_cds_mapping)


class CodingRuns:
    """A cached, compact transcript-coordinate index of genomic coding runs."""
    def __init__(self, cds):
        self.chrom = cds.iloc[0]['Seqid']
        self.strand = cds.iloc[0]['Strand']
        # Preserve exactly the established first-phase, tie-order and overlapping
        # CDS semantics. Build once per protein, never once per peptide.
        positions = np.asarray(_coding_nt_positions(cds), dtype=np.int64)
        self.length = len(positions)
        self.positions = None
        if self.length > 1 and np.any(np.diff(positions) * (1 if self.strand == '+' else -1) < 0):
            # Pathological overlapping/reordered CDSs: preserve arbitrary walks.
            self.positions = positions
        cuts = np.r_[0, np.flatnonzero(np.abs(np.diff(positions)) != 1) + 1, self.length]
        self.runs = [(int(a), int(b), int(positions[a]),
                      1 if self.strand == '+' else -1)
                     for a,b in zip(cuts[:-1],cuts[1:]) if a < b]
        self.ends = [run[1] for run in self.runs]

    def project(self, start, end):
        """Project 1-based inclusive amino-acid coordinates, preserving BED blocks."""
        if start < 1 or end < start or end*3 > self.length:
            return []
        left,right=(start-1)*3,end*3
        if self.positions is not None:
            selected=self.positions[left:right]
            cuts=np.r_[0,np.flatnonzero(np.abs(np.diff(selected)) != 1)+1,len(selected)]
            return sorted((int(selected[a:b].min()),int(selected[a:b].max()))
                          for a,b in zip(cuts[:-1],cuts[1:]) if a<b)
        regions=[]
        for index in range(bisect_right(self.ends,left),len(self.runs)):
            a,b,genomic,direction=self.runs[index]
            if a>=right:break
            x,y=max(a,left),min(b,right)
            first=genomic+direction*(x-a);last=genomic+direction*(y-a-1)
            regions.append((min(first,last),max(first,last)))
        return sorted(regions)


def attach_model_metadata(gtf_file, protein_fasta, mapping_file, output_dir):
    """Write a legacy-compatible metadata BED using cached CDS interval indexes.

    The downstream join retains only these records matching PepLocus loci, while
    preserving every additional genomic placement for ambiguity checks.
    """
    started=time.monotonic()
    print('Attaching protein metadata using cached CDS interval indexes...',flush=True)
    sequences={r.id:str(r.seq).replace('I','L').replace('*','') for r in SeqIO.parse(protein_fasta,'fasta')}
    annotation=gtf_to_df_with_genes(gtf_file,valid_protein_ids=sequences)
    mapping=mapping_file_to_df(mapping_file)
    mapping['Peptide']=mapping['Peptide'].str.replace('I','L')
    # Only Protein/Gene participates in this merge. Remove duplicate annotation
    # rows before joining rather than expanding every peptide over every CDS.
    peptides=get_gene_protein_specific_peps(annotation[['Protein','Gene']].drop_duplicates(),mapping)
    columns=['Peptide','Protein','Gene']+(['location'] if 'location' in peptides else [])
    peptides=peptides.drop_duplicates(subset=columns,keep='first')
    models=get_protein_to_cds_mapping(annotation,protein_ids=peptides['Protein'].dropna().unique())
    indexes={};records=[];seen=set()
    has_location='location' in peptides
    fields=['Peptide','Protein','Gene']+(['location'] if has_location else [])
    for values in peptides[fields].itertuples(index=False,name=None):
        peptide,protein,gene=values[:3];sequence=sequences.get(protein)
        if not sequence or protein not in models:continue
        index=indexes.get(protein)
        if index is None:
            index=indexes[protein]=CodingRuns(models[protein])
        for start in resolve_peptide_starts(sequence,peptide,values[3] if has_location else None):
            end=start+len(peptide)-1;regions=index.project(start,end)
            if not regions:continue
            key=(index.chrom,index.strand,peptide,protein,start,end,tuple(regions))
            if key in seen:continue
            seen.add(key)
            records.append(dict(seqid=index.chrom,strand=index.strand,peptide=peptide,
                                protein=protein,gene=gene,
                                type=classify_peptide_type(regions,start==1,end==len(sequence)),
                                pep_start=start,pep_end=end,prot_len=len(sequence),
                                previous_aa=sequence[start-2] if start>1 else '',genomic_regions=regions))
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    bed=output_dir/'peptides_mapped.bed';export_to_bed(records,bed)
    print(f'Attached {len(records)} protein metadata records in {time.monotonic()-started:.3f} s: {bed}',flush=True)
    return bed
