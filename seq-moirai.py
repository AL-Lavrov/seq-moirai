import modules.dna_rna_tools as drt

def run_dna_rna_tools(*seqs_and_action, seq_type=None):
    action = str(seqs_and_action[-1])
    seqs = seqs_and_action[:-1]

    if drt.is_nucleic_acid(seqs) is False:
        return None
    if seq_type is None:
        seq_type = drt.is_dna_or_rna(seqs)

    if action == 'reverse':
        processed_seq = [drt.reverse_seq(seq) for seq in seqs]
    elif action == 'transcribe' and seq_type:
        processed_seq = [drt.transcribe_seq(seq, seq_type) for seq in seqs]
    elif action == 'complement' and seq_type:
        processed_seq = [drt.complement_seq(seq, seq_type) for seq in seqs]
    elif action == 'reverse_complement' and seq_type:
        processed_seq = [drt.reverse_complement_seq(seq, seq_type) for seq in seqs]

    return processed_seq[0] if len(processed_seq) == 1 else processed_seq