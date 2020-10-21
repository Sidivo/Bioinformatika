from Bio import SeqIO
import numpy as np
from Bio.Data import CodonTable


# Funkcija nuskaityti faila
def read_fasta(filename):
    for sequence in SeqIO.parse(filename, "fasta"):
        return sequence
# Funkcija padaryti 6 frames ir tripletus
def get_triplets(sequence):
    frames = []
    frames.append([sequence.seq[i:i + 3] for i in range(0, len(sequence.seq), 3)])
    frames.append([sequence.seq[i:i + 3] for i in range(1, len(sequence.seq), 3)])
    frames.append([sequence.seq[i:i + 3] for i in range(2, len(sequence.seq), 3)])
    frames.append([sequence.seq.reverse_complement()[i:i + 3] for i in range(0, len(sequence.seq), 3)])
    frames.append([sequence.seq.reverse_complement()[i:i + 3] for i in range(1, len(sequence.seq), 3)])
    frames.append([sequence.seq.reverse_complement()[i:i + 3] for i in range(2, len(sequence.seq), 3)])
    return frames
# Funkcija rasti kodonus su ATG pradzia ir TAA arba TAG arba TGA pabaigomis
def find_codon(sequence):
    i = 0
    good_sequences = []
    while i < len(sequence):
        if sequence[i] == 'ATG':
            starting_place = i
            j = i
            while j < len(sequence):
                if sequence[j] == 'TAA' or sequence[j] == 'TAG' or sequence[j] == 'TGA':
                    ending_place = j
                    good_sequences.append(''.join(str(e) for e in sequence[starting_place:ending_place + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return good_sequences

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


