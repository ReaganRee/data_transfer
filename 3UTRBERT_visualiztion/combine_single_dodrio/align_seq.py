from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio.Seq import Seq


seq1 = Seq('UAGUCUCCCUCGGGGUCCGACAUACAGUGAUAGAGCAUGCUGUCCGCGCACCACAUGCAGCUCACCCGGCGGAUGCAAGUUCUCACGGAGUCGGGCGCGU')
seq2 = Seq('UCAGUGGUCCCGGGGUCCCCGGGACCUCUGUUGGCUGCGGCCCACUGCGGGCUGCAACCGCGGGCCGGGGCCGCGGGGAUGUGCAAAGGGCAGCGUCGGG')
alignments = pairwise2.align.globalxx(seq1, seq2)
for alignment in alignments:
    print(format_alignment(*alignment))


