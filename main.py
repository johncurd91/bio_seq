# Testing environment
from bio_seq import BioSeq
from utilities import read_text_file

rna_seq = BioSeq(seq_type='RNA', seq=read_text_file('test_file.txt'))

print(rna_seq.codon_usage('L'))
