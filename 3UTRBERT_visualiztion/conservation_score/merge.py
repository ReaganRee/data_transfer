from Bio import SeqIO

records_0_657 = list(SeqIO.parse(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/conservation_score/m6a_remapped_pos_0_657.fa",
        "fasta"))
records_658_1000 = list(SeqIO.parse(
        "/Users/reagan/Desktop/3UTRBERT_visualiztion/conservation_score/m6a_remapped_pos.fa",
        "fasta"))

records = records_0_657 + records_658_1000
SeqIO.write(records,"/Users/reagan/Desktop/3UTRBERT_visualiztion/conservation_score/m6a_remapped_0_1k.fa", 'fasta')
