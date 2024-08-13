# ORF-finder-using-python

## Description
This repository contains a Python script that identifies Open Reading Frames (ORFs) in a DNA sequence. The script analyzes all six reading frames and returns the ORFs along with their corresponding amino acid sequences.

## Features
- **Six Reading Frames Analysis**: The script checks all three forward and three reverse reading frames.
- **Codon Translation**: The script translates DNA codons to their respective amino acids.
- **Stop Codon Detection**: It detects ORFs by finding start (`ATG`) and stop codons (`TAA`, `TAG`, `TGA`).

## Usage
1. - **Clone the Repository**:
   ```bash
   git clone https://github.com/Kpmurshid/ORF-finder-using-python.git
   cd ORF-finder-using-python

2. - **Run the Script**:
   Use the following command to find ORFs in your DNA sequence:
   ```bash
     python find_orfs.py

3. - **Example**:
Here's how to use the script in Python:
   ```bash
      dna_sequence = "ATGAAATAGTGAATGCTAG"
      orf_dict = find_orfs(dna_sequence)
      
      for frame, orfs in orf_dict.items():
          print(f"ORF in Frame {frame}:")
          for idx, (orf, amino_acids) in enumerate(orfs, start=1):
              print(f"ORF {idx}: DNA: {orf}, Amino Acids: {amino_acids}")
          print()
