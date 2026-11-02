import importlib.util
from pathlib import Path
import random,csv,json

import pytest

from gapms.protein_mapper import ProteinMatcher


def mappings(path):
 with path.open() as handle:
  return {(r['peptide'],r['protein'],int(r['location'])) for r in csv.DictReader(handle,delimiter='\t')}


def oracle(peptides,proteins):
 expected=set()
 for peptide in peptides:
  found=False;q=peptide.replace('I','L')
  for name,sequence in proteins:
   text=sequence.upper().replace('I','L')
   for i in range(len(text)-len(q)+1):
    if text[i:i+len(q)]==q:expected.add((peptide,name,i+1));found=True
  if not found:expected.add((peptide,'UNMAPPED',0))
 return expected


@pytest.mark.parametrize('seed',range(10))
def test_randomized_exact_occurrence_oracle(tmp_path,seed):
 rng=random.Random(seed);proteins=[(f'p{i}',''.join(rng.choices('ACDEILKMP*X',k=150))) for i in range(20)]
 peptides={''.join(rng.choices('ACDEILKMP',k=rng.randint(1,15))) for _ in range(100)}
 peptides.update(['AAA','AAAA','IL','LL','I','L'])
 fasta=tmp_path/'proteins.fa';fasta.write_text(''.join(f'>{p}\n{s}\n' for p,s in proteins))
 matcher=ProteinMatcher(peptides);out=tmp_path/'mapping.tsv';matcher.map(fasta,out)
 assert mappings(out)==oracle(peptides,proteins)


def test_automaton_reuse_and_boundaries(tmp_path):
 matcher=ProteinMatcher(['AA','AIL','ALL','MAK','QQQ']);automaton=matcher.automaton
 for i,sequence in enumerate(['AAAA*AILXALL','ALL']):
  fasta=tmp_path/f'{i}.fa';fasta.write_text('>one\n'+sequence+'\n>two\nMA*AK\n')
  out=tmp_path/f'{i}.tsv';matcher.map(fasta,out)
  assert matcher.automaton is automaton
  assert mappings(out)==oracle(matcher.peptides,[('one',sequence),('two','MA*AK')])
  summary=json.loads(out.with_suffix('.summary.json').read_text())
  assert summary['automaton_builds']==1 and summary['automaton_searches']==i+1


def test_fragpipe_normalization_and_no_header_peptide(tmp_path):
 peptides=tmp_path/'peptides.tsv';peptides.write_text('Peptide Sequence\nK.AI[Oxidation]LK.R\nAILK\n')
 matcher=ProteinMatcher.from_file(peptides)
 assert matcher.peptides==('AILK',)
 fasta=tmp_path/'p.fa';fasta.write_text('>p description\nALLK\n')
 out=tmp_path/'mapping.tsv';matcher.map(fasta,out)
 assert mappings(out)=={('AILK','p',1)}
 with pytest.raises(ValueError,match='already exists'):matcher.map(fasta,out)


def test_duplicate_identifiers_and_empty_input(tmp_path):
 matcher=ProteinMatcher(['AAA']);fasta=tmp_path/'p.fa';out=tmp_path/'out.tsv'
 fasta.write_text('>p\nAAA\n>p\nAAA\n')
 with pytest.raises(ValueError,match='duplicate'):matcher.map(fasta,out)
 assert not out.exists()
 fasta.write_text('')
 with pytest.raises(ValueError,match='empty'):matcher.map(fasta,out)
