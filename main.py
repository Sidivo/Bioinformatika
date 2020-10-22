from Bio import SeqIO
import numpy as np
from Bio.Data import CodonTable


# Funkcija nuskaityti faila
def read_fasta(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        return seq_record
# Funkcija padaryti 6 frames ir tripletus
def get_triplets(seq_record):
    frames = []
    frames.append([seq_record.seq[i:i + 3] for i in range(0, len(seq_record.seq), 3)])
    frames.append([seq_record.seq[i:i + 3] for i in range(1, len(seq_record.seq), 3)])
    frames.append([seq_record.seq[i:i + 3] for i in range(2, len(seq_record.seq), 3)])
    frames.append([seq_record.seq.reverse_complement()[i:i + 3] for i in range(0, len(seq_record.seq), 3)])
    frames.append([seq_record.seq.reverse_complement()[i:i + 3] for i in range(1, len(seq_record.seq), 3)])
    frames.append([seq_record.seq.reverse_complement()[i:i + 3] for i in range(2, len(seq_record.seq), 3)])
    return frames
# Funkcija rasti kodonus su ATG pradzia ir TAA arba TAG arba TGA pabaigomis
def find_codon(seq_record):
    i = 0
    codon_list = []
    while i < len(seq_record):
        if seq_record[i] == 'ATG':
            start_pos = i
            j = i
            while j < len(seq_record):
                if seq_record[j] == 'TAA' or seq_record[j] == 'TAG' or seq_record[j] == 'TGA':
                    end_pos = j
                    codon_list.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codon_list

def find_codons(triplets):
    codons = []
    for triplet in triplets:
        codons.append(find_codon(triplet))
    return codons

def longest_codon(codons):
    return max(codons, key=len)


if __name__ == '__main__':

    #nuskaitom
    record = read_fasta('bacterial1.fasta')
    #gaunam tripletus
    triplets = get_triplets(record)
    #random kodonus
    for triplet in triplets:
        len(triplet)
    codons = find_codons(triplets)
    #sujungiam i viena
    codons = np.concatenate(codons)
    #randam didziausia
    max_codon = longest_codon(codons)
    #atrusiuojam didziausius, kurie ilgesni nei 100
    longest_codons = []
    for codon in codons:
        if len(codon) >= 100:
            longest_codons.append(codon)

    print(longest_codons)


