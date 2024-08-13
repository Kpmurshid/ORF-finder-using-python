from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class ORFFinder:
    """A class to identify Open Reading Frames (ORFs) in DNA sequences."""
    
    def __init__(self, sequence_record, sequence_type='DNA'):
        self.sequence_record = sequence_record
        self.sequence_type = sequence_type
        self.start_codon = 'ATG' if sequence_type == 'DNA' else 'AUG'
        self.orfs = {'plus': [], 'minus': []}
    
    def _find_orfs_in_strand(self, genetic_code_table, is_reverse, offset, min_length):
        """Find ORFs in a single strand of the sequence."""
        orf_records = []
        strand_seq = self.sequence_record.seq.reverse_complement() if is_reverse else self.sequence_record.seq
        start_positions = self._find_start_positions(strand_seq)
        
        for index, start_pos in enumerate(start_positions, start=1):
            protein_seq = strand_seq[start_pos:].translate(table=genetic_code_table, to_stop=True)
            if len(protein_seq) >= min_length:
                orf_id = f"{self.sequence_record.id}_{index + offset}"
                end_pos = start_pos + len(protein_seq) * 3
                description = self._generate_description(start_pos, end_pos, is_reverse)
                
                orf_record = SeqRecord(protein_seq, id=orf_id, description=description)
                orf_records.append(orf_record)
        
        return orf_records
    
    def _find_start_positions(self, sequence):
        """Identify all start positions of the ORFs."""
        positions = []
        pos = 0
        while True:
            pos = sequence.find(self.start_codon, pos)
            if pos == -1:
                break
            positions.append(pos)
            pos += len(self.start_codon)
        return positions
    
    def _generate_description(self, start, end, is_reverse):
        """Generate a description for an ORF based on its position and orientation."""
        length = len(self.sequence_record.seq)
        if is_reverse:
            start = length - end + 1
            end = length - start + len(self.start_codon)
            return f"[{start} - {end}] (REVERSE SENSE) {self.sequence_record.description}"
        else:
            return f"[{start + 1} - {end}] {self.sequence_record.description}"
    
    def find_orfs(self, genetic_code_table, min_length):
        """Find ORFs in both strands of the sequence."""
        self.orfs['plus'] = self._find_orfs_in_strand(genetic_code_table, is_reverse=False, offset=0, min_length=min_length)
        self.orfs['minus'] = self._find_orfs_in_strand(genetic_code_table, is_reverse=True, offset=len(self.orfs['plus']), min_length=min_length)
        return self.orfs

# Example usage
if __name__ == "__main__":
    # Example DNA sequence
    dna_seq = Seq("ATGCGTACATCGCAGATGCAGTACGAGGACTAGCATCACA")
    record = SeqRecord(dna_seq, id="example_sequence", description="Example DNA sequence")

    # Create ORF Finder instance
    orf_finder = ORFFinder(record, sequence_type='DNA')

    # Find ORFs with genetic code table 1 and minimum length of 10 amino acids
    orf_results = orf_finder.find_orfs(genetic_code_table=1, min_length=10)

    # Print ORF results
    for strand, orfs in orf_results.items():
        for orf in orfs:
            print(f"{strand.capitalize()} Strand ORF:", orf.id, orf.description, orf.seq)

    print("Total ORFs found:", len(orf_results['plus']) + len(orf_results['minus']))
