import bases


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
