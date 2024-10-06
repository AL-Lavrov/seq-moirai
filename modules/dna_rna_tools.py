from . import bases
from typing import Tuple, Union, TypedDict


class Options(TypedDict):
    seq_type: Union[None, str]
    frame: int


def get_data(*seqs_and_action: str,
             **options: Options
             ) -> Tuple[Tuple[str], str, Options]:

    action = str(seqs_and_action[-1])
    options = add_defaults(**options)
    seqs = seqs_and_action[:-1]
    if is_nucleic_acid(seqs) is False:
        raise ValueError("Not a nucleic acid")
    if options["seq_type"] is None:
        options["seq_type"] = is_dna_or_rna(seqs)

    return seqs, action, options


def add_defaults(**options: Options) -> Options:
    default_options = {"seq_type": None, "frame": 1}
    default_options.update(**options)
    return default_options


def is_nucleic_acid(seqs: tuple) -> bool:
    seqs = tuple(seq.upper() for seq in seqs)
    return all(seq and set(seq).issubset(bases.nucleotides) for seq in seqs)


def is_dna_or_rna(seqs: tuple) -> str:
    seqs = tuple(seq.upper() for seq in seqs)
    if all(set(seq).issubset(bases.dna_bases) for seq in seqs):
        return "dna"
    elif all(set(seq).issubset(bases.rna_bases) for seq in seqs):
        return "rna"
    else:
        raise ValueError("Mixed input of DNA and RNA")


def transcribe(seq: str, **options: Options) -> str:
    seq_type = options["seq_type"]
    if seq_type == "dna":
        return seq.replace("t", "u").replace("T", "U")
    if seq_type == "rna":
        return seq.replace("u", "t").replace("U", "T")


def reverse(seq: str, **options: Options) -> str:
    return seq[::-1]


def complement(seq: str, **options: Options) -> str:
    seq_type = options["seq_type"]
    if seq_type == "rna":
        return seq.translate(str.maketrans(bases.rna_base_pairs))
    if seq_type == "dna":
        return seq.translate(str.maketrans(bases.dna_base_pairs))


def reverse_complement(seq: str, **options: Options) -> str:
    return reverse(complement(seq, **options), **options)


def count_gc(seq: str, **options: Options) -> str:
    seq = seq.upper()
    gc_content = (seq.count("G") + seq.count("C")) / len(seq)
    return str(gc_content)


def translate(seq: str, **options: Options) -> str:
    protein, frame = "", options["frame"]
    seq = seq.upper()[frame - 1:]
    if options["seq_type"] == "dna":
        seq = transcribe(seq, **options)
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) == 3:
            protein += bases.codons[codon]
    return protein
