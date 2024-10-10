from . import bases
from typing import Tuple, Union, TypedDict


class Options(TypedDict):
    '''
    Default structure of dna_rna_tools options.

    seq_type (str): Sequence type of input nucleic acids. 
        'rna' and 'dna' are accepted by downstream functions.
        None value triggers get_data function to check nucleic
        acid type in input.
    frame (int): Number of reading frame to translate. Accepted 
        values are 1, 2, 3.
    '''
    seq_type: Union[None, str]
    frame: int


def get_data(*seqs_and_action: str,
             **options: Options
             ) -> Tuple[Tuple[str], str, Options]:
    '''
    Orginizes data to fit run_dna_rna_tools function structure.

    Splits positional arguments into sequences and action, validates
    sequences and forms a complete class Options dictionary.

    Arguments
    ---------
    *seqs_and_action (str): any number of nucleic, followed by action.

    **options: Additional keyword arguments.
    - seq_type (str): Type of nucleic acid to be processed by
    rna_dna_tools.
    - frame (int): Specifies the reading frame 
    for the 'translate' action.
    Returns
    -------

    Tuple[Tuple[str], str, Options]: structured data which can be processed 
    by run_rna_dna_tools.
    '''
    action = str(seqs_and_action[-1])
    options = add_defaults(**options)
    seqs = seqs_and_action[:-1]
    if is_nucleic_acid(seqs) is False:
        raise ValueError("Not a nucleic acid")
    if options["seq_type"] is None:
        options["seq_type"] = is_dna_or_rna(seqs)

    return seqs, action, options


def add_defaults(**options: Options) -> Options:
    '''
    Add default values to class Options dictionary if some 
    are not specified.

    Arguments
    ---------
    **options: Additional keyword arguments.

    - seq_type (str): a sequence type specified.
    - frame (int): a reading frame specified.

    Returns
    -------
    dict (Options): A dictionary class Options with user specified
    values and defualt values if one wasn't provided. Default values:
    - seq_type: None
    - frame: 1
    '''
    default_options = {"seq_type": None, "frame": 1}
    default_options.update(**options)
    return default_options


def is_nucleic_acid(seqs: tuple) -> bool:
    '''
    Checks if all of sequences passed are valid nucleic acids
     
    Valid nucleic acids contain only {A, T, G, C, U} bases, 
    do not contain T and U simultaneously, must not be
    blank strings.

    Arguments
    ---------
    seqs (Tuple[str]): one or more sequences.

    Returns
    --------
    bool: True if condtions are satisfied for all sequences.
    False otherwise.
    '''
    seqs = tuple(seq.upper() for seq in seqs)
    return all(seq and set(seq).issubset(bases.nucleotides) for seq in seqs)


def is_dna_or_rna(seqs: tuple) -> str:
    '''
    Checks the type of nucleic acids sequences passed.
     
    Arguments
    ---------
    seqs (Tuple[str]): one or more sequences.

    Returns
    --------
    str
    '''
    seqs = tuple(seq.upper() for seq in seqs)
    if all(set(seq).issubset(bases.dna_bases) for seq in seqs):
        return "dna"
    elif all(set(seq).issubset(bases.rna_bases) for seq in seqs):
        return "rna"
    else:
        raise ValueError("Mixed input of DNA and RNA")


def transcribe(seq: str, **options: Options) -> str:
    '''
    Converts nucleic acid sequence in sequence of opposite type.

    Arguments
    ---------
    seq (str): The nucleic acid sequence which must be transcribed.

    **options (Options): Keyword arguments.

    - seq_type (str): Specifies direction of transcription. 
    - Must be 'dna' or 'rna'.
    
    Returns
    -------
    str: nucleic acid sequence of opposite type.
    '''
    seq_type = options["seq_type"]
    if seq_type == "dna":
        return seq.replace("t", "u").replace("T", "U")
    if seq_type == "rna":
        return seq.replace("u", "t").replace("U", "T")


def reverse(seq: str, **options: Options) -> str:
    '''
    Reverses a sequence of nucleic acid.

    Arguments
    ---------
    seq (str): a sequence of nucleic acid to reverse.

    **options (Options): additional keyword arguments.

    - Can be passed, but won't be used

    Returns
    -------
    str: passed sequnce in reverse direction.
    '''
    return seq[::-1]


def complement(seq: str, **options: Options) -> str:
    '''
    Generates 3'-5' complementary sequence.

    DNA complement is generated for input DNA sequences, and RNA
    complement is generated for RNA sequences.

    Arguments
    ---------
    seq (str): The nucleic acid sequence for which complement 
    must be generated.

    **options (Options): Keyword arguments.

    - seq_type (str): specifies a complementary rules applied for
    input sequence ('dna' or 'rna').
    
    Returns
    -------
    str: complementary strand for input sequence in 3'-5' direction.
    '''
    seq_type = options["seq_type"]
    if seq_type == "rna":
        return seq.translate(str.maketrans(bases.rna_base_pairs))
    if seq_type == "dna":
        return seq.translate(str.maketrans(bases.dna_base_pairs))


def reverse_complement(seq: str, **options: Options) -> str:
    '''
    Generates 5'-3' complementary sequence.

    DNA complement is generated for input DNA sequences, and RNA
    complement is generated for RNA sequences. This function combines 
    the operations of complementing and reversing.

    Arguments
    ---------
    seq (str): The nucleic acid sequence for which complement 
    must be generated.

    **options (Options): Keyword arguments.

    - seq_type (str): specifies a complementary rules applied for
    - input sequence. Must be 'dna' or 'rna'.
    
    Returns
    -------
    str: complementary strand for input sequence in 5'-3' direction.
    '''
    return reverse(complement(seq, **options), **options)


def count_gc(seq: str, **options: Options) -> str:
    '''
    Counts fraction of G and C bases in input sequence.

    Arguments
    ---------
    seq (str): a sequence for which GC content must be count.

    **options (Options): Additional keyword arguments.

    - Can be passed, but won't be used

    Returns
    -------
    (str) - a fraction of G and C in input sequence.
    
    Notes
    -----
    Computed GC-content is converted to string to unify output
    of all module's functions.
    '''
    seq = seq.upper()
    gc_content = (seq.count("G") + seq.count("C")) / len(seq)
    return str(gc_content)


def translate(seq: str, **options: Options) -> str:
    '''
    Translates nucleic acid sequence into the protein it encodes.

    Translation occurs according to the standard genetic code in the
    specified reading frame. The translation starts at the first, 
    second, or third base, depending on the chosen reading frame, and 
    any bases that do not form a complete triplet are truncated.

    Arguments
    ---------
    seq (str): A nucleic acid sequence to be translated to protein.

    **options: Keyword arguments.

    - seq_type (str): Type of sequence to be translated ('rna' or 'dna')
    - frame (int): A reading frame in which translation will occur (1, 2, 3).

    Returns
    -------
    str: Translated protein sequence in one-letter notation. 
    '''
    protein, frame = "", options["frame"]
    seq = seq.upper()[frame - 1:]
    if options["seq_type"] == "dna":
        seq = transcribe(seq, **options)
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) == 3:
            protein += bases.codons[codon]
    return protein
