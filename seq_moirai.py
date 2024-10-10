import modules.dna_rna_tools as drt
import modules.filter_fastq as ffq
from typing import Union, Tuple, Dict


def run_dna_rna_tools(*seqs_and_action: str,
                      **options: drt.Options
                      ) -> Union[str, list]:
    '''
    Transforms nucleic acid sequences.

    This function validates the input (nucleic acid sequences and
    the action to be performed) and applies the corresponding
    transformation to each sequence.

    Arguments
    ---------

    *seqs_and_action (str): Any number of nucleic acid sequences
        followed by the desired action as the final argument.
        All sequences must be of the same type (either RNA or DNA),
        and degenerate bases (Y, N, etc.) are not allowed.
        Accepted actions:
        - 'transcribe' - Converts DNA to RNA.
        - 'reverse' - Reverses the sequence (3' to 5' direction).
        - 'complement' - Generates the complementary sequence
          (in the 3' to 5' direction).
        - 'reverse-complement' - Returns the reverse-complementary
          sequence (in the 5' to 3' direction).
        - 'count_gc' - Returns a string with the GC content of
          the sequence.
        - 'translate' - Translates the sequence into a protein.
          Start codons are marked as 'M!', and stop codons as '*'.
          The reading frame is determined by the `frame` argument.
          Bases not forming a complete triplet are truncated.

    **options: Additional keyword arguments.
    - seq_type (str): Typically not required, but use 'rna' for RNA
    nputs lacking uridine (U) bases. This overrides mixed input control.
    - frame (int): Can be 1, 2, or 3. Specifies the reading frame
    for the 'translate' action. Default is 1.

    Returns
    -------
    A string with the transformed sequence if one sequence is passed,
    or a list of strings if two or more sequences are passed.
    '''
    seqs, action, options = drt.get_data(*seqs_and_action, **options)
    action_map = {"transcribe": drt.transcribe,
                  "reverse": drt.reverse,
                  "reverse_complement": drt.reverse_complement,
                  "complement": drt.complement,
                  "count_gc": drt.count_gc,
                  "translate": drt.translate}
    processed_seq = [action_map[action](seq, **options) for seq in seqs]
    return processed_seq[0] if len(processed_seq) == 1 else processed_seq


def filter_fastq(fastq_entry: Dict[str, Tuple[str]],
                 gc_bounds: Union[float, Tuple[float]] = (0, 100),
                 length_bounds: Union[float, Tuple[float]] = (0, 2**32),
                 quality_threshold: int = 0,
                 ) -> Dict[str, Tuple[str]]:
    '''
    Filters a FASTQ entry based on specified GC-content, length
    and quality bounds.

    The function selects sequences that lie within the specified bounds
    of GC-content, length, and have a quality score higher than the
    specified threshold, and returns a new dictionary
    consisting of the sequences that pass the filtration.

    Arguments
    ---------

    fastq_entry (Dict[str, Tuple[str]]): Fastq entry which is represented as
      dictionary with structure {'name' : ('sequence', 'quality')}.
      - 'name' - sequence ID
      - 'sequence' - read sequence consisting of {A, T, G, C}.
      - 'quality' - A string with ASCII-coded quality of each base read.

    gc_bounds Union[float, Tuple[float]]: Bounds of accepted GC-content.
    Can be given as single number (treated as upper limit) or as an
    interval (tuple). Defualt is (0, 100).

    length_bounds Union[float, Tuple[float]]: Bounds of accepted length.
    Can be given as single number (treated as upper limit) or as an
    interval (tuple). Default is (0, 2**32).

    quality_threshold (int): The lower limit of Q-score. Default is 0.

    Returns
    -------
    (Dict[str, Tuple[str]]): Dictionary of the same structure as input
    dictionary, but containing only sequences which passed filtration.
    '''
    processed_fastq = {}
    for name, seq_info in fastq_entry.items():
        seq, seq_quality = seq_info
        if (
            ffq.filter_gc(seq, gc_bounds)
            and ffq.filter_length(seq, length_bounds)
            and ffq.filter_quality(seq_quality, quality_threshold)
          ):
            processed_fastq.update({name: seq_info})
    return processed_fastq
