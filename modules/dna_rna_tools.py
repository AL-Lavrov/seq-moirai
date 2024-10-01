import bases


def run_dna_rna_tools(*seqs_and_action, seq_type=None):
    action = str(seqs_and_action[-1])
    seqs = seqs_and_action[:-1]

    if is_nucleic_acid(seqs) is False:
        return None
    if seq_type is None:
        seq_type = is_dna_or_rna(seqs)

    if action == 'reverse':
        processed_seq = [reverse_seq(seq) for seq in seqs]
    elif action == 'transcribe' and seq_type:
        processed_seq = [transcribe_seq(seq, seq_type) for seq in seqs]
    elif action == 'complement' and seq_type:
        processed_seq = [complement_seq(seq, seq_type) for seq in seqs]
    elif action == 'reverse_complement' and seq_type:
        processed_seq = [reverse_complement_seq(seq, seq_type) for seq in seqs]

    return processed_seq[0] if len(processed_seq) == 1 else processed_seq


def is_nucleic_acid(seqs: tuple) -> bool:
    seqs = tuple(seq.upper() for seq in seqs)
    return all(seq and set(seq).issubset(bases.nucleotides) for seq in seqs)


def is_dna_or_rna(seqs: tuple) -> str:
    seqs = tuple(seq.upper() for seq in seqs)
    if all(set(seq).issubset(bases.dna_bases) for seq in seqs):
        return 'dna'
    elif all(set(seq).issubset(bases.rna_bases) for seq in seqs):
        return 'rna'
    else:
        return None


def transcribe_seq(seq: str, seq_type: str) -> str:
    if seq_type == 'dna':
        return seq.replace('t', 'u').replace('T', 'U')
    if seq_type == 'rna':
        return seq.replace('u', 't').replace('U', 'T')


def reverse_seq(seq: str) -> str:
    return seq[::-1]


def complement_seq(seq: str, seq_type: str) -> str:
    if seq_type == "rna":
        return seq.translate(str.maketrans(bases.rna_base_pairs))
    if seq_type == "dna":
        return seq.translate(str.maketrans(bases.dna_base_pairs))


def reverse_complement_seq(seq: str, seq_type: str) -> str:
    return reverse_seq(complement_seq(seq, seq_type))
