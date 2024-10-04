import bases

def orginize_data(*seqs_and_action, **options):
    print(options)
    seqs = seqs_and_action[:-1]
    action = str(seqs_and_action[-1])
    if is_nucleic_acid(seqs) is False:
        raise ValueError('Not a nucleic acid')
    if 'seq_type' not in options:
        options['seq_type'] = is_dna_or_rna(seqs)
    return seqs, action, {'seq_type': options['seq_type'],
        'frame': options['frame']}


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


def reverse_seq(seq, **kwargs) -> str:
    return seq[::-1]


def complement_seq(seq, **options) -> str:
    seq_type = options['seq_type']
    if seq_type == "rna":
        return seq.translate(str.maketrans(bases.rna_base_pairs))
    if seq_type == "dna":
        return seq.translate(str.maketrans(bases.dna_base_pairs))


def reverse_complement_seq(seq: str, seq_type: str) -> str:
    return reverse_seq(complement_seq(seq, seq_type))

action_map = {'transcribe': transcribe_seq,
              'reverse': reverse_seq,
              'reverse_complement': reverse_complement_seq,
              'complement': complement_seq}

print(orginize_data('GAATAGGGgc', 'uuuaaaa', 'complement'))
