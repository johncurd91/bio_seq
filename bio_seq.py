from bio_structs import DNA_codons, RNA_codons, NUCLEOTIDE_BASE
from collections import Counter
import random


class BioSeq:
    """DNA sequence class; seq='ATCG', seq_type='DNA', label='No Label' by default"""

    def __init__(self, seq='ATCG', seq_type='DNA', label='No Label'):
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.label = label
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not appear to be a valid {self.seq_type} string."

    # DNA Toolkit functions:
    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n" \
               f"[Sequence]: {self.seq}\n" \
               f"[Biotype]: {self.seq_type}\n" \
               f"[Length]: {len(self.seq)}"

    def gen_rand_seq(self, length=10, seq_type='DNA'):
        """Generate a random DNA sequence; length=10 by default"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, 'Randomly generated sequence')

    def nucleotide_freq(self):
        """Count nucleotides in given sequence. Return a dictionary"""
        return dict(Counter(self.seq))

    def transcription(self):
        """DNA > RNA transcription, replacing T with U"""
        return self.seq.replace('T', 'U')

    def reverse_transcription(self):
        """RNA > DNA reverse transcription, replacing U with T"""
        if self.seq_type == 'DNA':
            return self.seq.replace('U', 'T')
        return 'Not a DNA sequence'

    def complement(self):
        """Return the complement of a given nucleotide sequence."""
        if self.seq_type == 'DNA':
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)

    def reverse_complement(self):
        """Return the REVERSE complement of a given nucleotide sequence."""
        if self.seq_type == 'DNA':
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """Calculates GC content in DNA/RNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)

    def gc_content_subseq(self, k=20):
        """GC content in a DNA/RNA sub-sequence; length k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        """Translate nucleotide sequence into amino acid sequence"""
        if self.seq_type == 'DNA':
            codon_table = DNA_codons
        else:
            codon_table = RNA_codons

        protein = [codon_table[self.seq[pos:pos + 3]]
                   for pos in range(init_pos, len(self.seq) - 2, 3)]
        return ''.join(protein)

    def codon_usage(self, amino_acid):
        """Provide frequency of each codon encoding a given amino acid in a nucleotide seq."""
        temp_lst = []
        if self.seq_type == 'DNA':
            codon_table = DNA_codons
        else:
            codon_table = RNA_codons

        for i in range(0, len(self.seq) - 2, 3):
            if codon_table[self.seq[i:i + 3]] == amino_acid:
                temp_lst.append(self.seq[i:i + 3])

        freq_dict = dict(Counter(temp_lst))
        total_wight = sum(freq_dict.values())
        for seq in freq_dict:
            freq_dict[seq] = round(freq_dict[seq] / total_wight, 2)

        return freq_dict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including rev. complement."""
        tmp_seq = BioSeq(self.reverse_complement(), self.seq_type)
        frames = [self.translate_seq(0),
                  self.translate_seq(1),
                  self.translate_seq(2),
                  tmp_seq.translate_seq(0),
                  tmp_seq.translate_seq(1),
                  tmp_seq.translate_seq(2)]
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an amino acid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating amino acids if '_' STOP codon was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if 'M' START codon was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, start_pos=0, end_pos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        if end_pos > start_pos:
            tmp_seq = BioSeq(
                self.seq[start_pos:end_pos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
            del tmp_seq
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

    def inverted_tandem_repeat(self, itr_len):
        """ Searches for inverted tandem repeats of given length. """
        for i in range(len(self.seq) - itr_len - 1):
            sub_seq = self.seq[i:i + itr_len]
            if sub_seq == sub_seq[::-1]:
                if self.seq_type == 'DNA':
                    mapping = str.maketrans('ATCG', 'TAGC')
                else:
                    mapping = str.maketrans('AUCG', 'UAGC')
                sub_seq_rev = sub_seq.translate(mapping)[::-1]

                for j in range(len(self.seq) - itr_len):
                    sub_seq_2 = self.seq[j:j + itr_len]
                    if sub_seq_2 == sub_seq_rev:
                        if j > i+8:
                            print(f"{sub_seq} at pos {i}:{i+8} = {sub_seq_2} at pos {j}:{j+8}.")



