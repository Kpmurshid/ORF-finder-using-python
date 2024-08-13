#Input sequence
dna_sequence = "ATGAAATAGTGAATGCTAG"

#Function for finding ORF
def find_orfs(dna_sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = {frame: [] for frame in range(1, 7)}
    codon_to_aa = {
        "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
        "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
        "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
        "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
    }

    for frame in range(3):
        frame_sequence = dna_sequence[frame:]
        for i in range(0, len(frame_sequence) - 2, 3):
            codon = frame_sequence[i:i+3]
            if codon == start_codon:
                amino_acid_sequence = ""
                for j in range(i, len(frame_sequence), 3):
                    codon = frame_sequence[j:j+3]
                    if codon in stop_codons:
                        orfs[frame + 1].append((frame_sequence[i:j+3], amino_acid_sequence))
                        break
                    amino_acid = codon_to_aa.get(codon, "?")
                    amino_acid_sequence += amino_acid

    reverse_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_sequence = "".join([reverse_complement[base] for base in dna_sequence[::-1]])

    for frame in range(3):
        frame_sequence = reverse_sequence[frame:]
        for i in range(0, len(frame_sequence) - 2, 3):
            codon = frame_sequence[i:i+3]
            if codon == start_codon:
                amino_acid_sequence = ""
                for j in range(i, len(frame_sequence), 3):
                    codon = frame_sequence[j:j+3]
                    if codon in stop_codons:
                        orfs[frame + 4].append((reverse_sequence[i:j+3], amino_acid_sequence))
                        break
                    amino_acid = codon_to_aa.get(codon, "?")
                    amino_acid_sequence += amino_acid

    return orfs

#Finding ORF
orf_dict = find_orfs(dna_sequence)

for frame, orfs in orf_dict.items():
    print(f"ORF in Frame {frame}:")
    for idx, (orf, amino_acids) in enumerate(orfs, start=1):
        print(f"ORF {idx}: DNA: {orf}, Amino Acids: {amino_acids}")
    print()


