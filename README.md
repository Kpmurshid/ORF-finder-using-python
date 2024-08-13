# ORF Finder

## Description

The ORF Finder is a Python script designed to detect Open Reading Frames (ORFs) in DNA sequences. This script uses the Biopython library to identify ORFs in both the plus and minus strands of a DNA sequence, based on user-defined parameters. The results include the position and translation of each ORF, which can be useful for genomic studies and protein prediction.

## Features

- Detect ORFs on both plus and minus DNA strands.
- Specify genetic code tables for translation.
- Filter ORFs by minimum length.
- Outputs ORF details including positions and translations.

## Requirements

- Python 3.x
- Biopython library

## Installation

1. Install the Biopython library:

   ```bash
   pip install biopython
