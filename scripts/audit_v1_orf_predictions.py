"""Strict protein identity, with proven deterministic IUPAC exceptions only."""
from pathlib import Path
import csv,json,re
from itertools import product
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values
import pysam

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('jobs_manifest',type=Path)
args=parser.parse_args()
root=args.jobs_manifest.resolve().parent;records=json.loads(args.jobs_manifest.read_text());report=[]
for record in records:
 out=Path(record['output'])
 original={r.id:str(r.seq).upper().rstrip('*') for r in SeqIO.parse(out/'transcripts/assembled_transcripts.fa.transdecoder.pep','fasta')}
 with (out/'orfs/orf_peptide_evidence.tsv').open() as stream:model_identity={r['model_id']:r['original_transcript_id'] for r in csv.DictReader(stream,delimiter='\t')}
 total=0;exceptions=[]
 for r in SeqIO.parse(out/'orfs/candidate_orfs.faa','fasta'):
  total+=1;tid=model_identity[r.id];predicted=original[tid];actual=str(r.seq)
  if predicted!=actual:exceptions.append(dict(model=r.id,transdecoder_id=tid,predicted=predicted,actual=actual))
 needed={r['transdecoder_id'] for r in exceptions};cds={r.id:str(r.seq).upper() for r in SeqIO.parse(out/'transcripts/assembled_transcripts.fa.transdecoder.cds','fasta') if r.id in needed}
 blocks={tid:[] for tid in needed}
 with (out/'orfs/candidate_orfs.gtf').open() as stream:
  for line in stream:
   if line.startswith('#') or not line.strip():continue
   fields=line.rstrip('\n').split('\t')
   if fields[2]!='CDS':continue
   tid=re.search(r'transcript_id "([^\"]+)"',fields[8]).group(1)
   if tid in blocks:blocks[tid].append((fields[0],int(fields[3])-1,int(fields[4]),fields[6],fields[7]))
 proven=[]
 with pysam.FastaFile(record['assembly'],filepath_index=str(next((out/'work/reference').glob('*.fai')))) as genome:
  for item in exceptions:
   tid=item['transdecoder_id'];segments=sorted(blocks[tid],key=lambda s:s[1]);strand=segments[0][3]
   nt=''.join(genome.fetch(chrom,a,b) for chrom,a,b,_,phase in segments).upper()
   if strand=='-':nt=str(Seq(nt).reverse_complement())
   phase=segments[0 if strand=='+' else -1][4];nt=nt[int(phase) if phase!='.' else 0:]
   assert nt==cds[tid],'transcript/genome nucleotide paths disagree'
   translated=str(Seq(nt[:len(nt)//3*3]).translate()).rstrip('*')
   assert translated==item['actual'] and len(item['predicted'])==len(item['actual'])
   resolved=[]
   for index,(predicted,actual) in enumerate(zip(item['predicted'],item['actual'])):
    if predicted==actual:continue
    codon=nt[index*3:index*3+3]
    assert predicted=='X' and set(codon)-set('ACGT'),'not an IUPAC masking difference'
    alternatives={str(Seq(''.join(bases)).translate()) for bases in product(*(ambiguous_dna_values[base] for base in codon))}
    assert alternatives=={actual},(codon,alternatives,actual)
    resolved.append(dict(amino_acid_position=index,transdecoder_aa=predicted,genome_aa=actual,codon=codon,possible_amino_acids=sorted(alternatives)))
   proven.append(dict(model=item['model'],transdecoder_id=tid,nucleotide_path_identical=True,deterministic_iupac_resolutions=resolved))
 row=dict(crop=record['crop'],status='validated',orfs=total,literal_protein_matches=total-len(exceptions),proven_iupac_masking_cases=len(proven),unexplained_disagreements=0,exceptions=proven)
 report.append(row);print(row['crop'],total,len(proven),'proven IUPAC cases',flush=True)
 (root/'orf_translation_audit_resolved.json').write_text(json.dumps(report,indent=2)+'\n')
assert len(report)==len(records) and all(row['status']=='validated' and row['unexplained_disagreements']==0 for row in report)
