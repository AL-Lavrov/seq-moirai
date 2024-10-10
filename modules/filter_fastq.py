from typing import Union, Tuple


def filter_gc(seq: str, gc_bounds: Union[float, Tuple[float]]) -> bool:
    '''
    Checks if nucleic acid sequence fits the specified GC-content.

    Arguments
    ---------
    seq (str): A nucleic acid sequence to be analyzed.

    gc_bounds Union[float, Tuple[float]]: Bounds of accepted GC-content.
    Can be given as single number (treated as upper limit) or as an
    interval (tuple).

    Returns
    -------
    bool: True if GC-content of sequence lies within the bounds. False otherwise.
    '''
    if isinstance(gc_bounds, Union[float, int]):
        gc_bounds = (0, gc_bounds)
    seq = seq.upper()
    gc_content = (seq.count("G") + seq.count("C")) / len(seq)
    return gc_bounds[0] <= gc_content * 100 <= gc_bounds[1]


def filter_length(seq: str, length_bounds: Union[int, Tuple[int]]) -> bool:
    '''
    Checks if nucleic acid sequence fits the specified lenght bounds.

    Arguments
    ---------
    seq (str): A nucleic acid sequence to be analized.

    length_bounds Union[float, Tuple[float]]: Bounds of accepted length.
    Can be given as single number (treated as upper limit) or as an
    interval (tuple).

    Returns
    -------
    bool: True if length of sequence lies within the bounds. False otherwise.
    '''
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def filter_quality(seq_quality: str, quality_threshold: int) -> bool:
    '''
    Checks if mean quality of sequence fits specified threshold of Q-score.

    Converts ASCII-coded quality of each base of sequence to Phred33 score, 
    evaluates the mean quality, and compares it with specified treshold.  

    Arguments
    ---------
    seq (str): A string with ASCII-coded quality of each base read.

    quality_threshold (int): The lower limit of Q-score.

    Returns
    -------
    bool: True if mean quality of sequence is greater or equal to the
    threshold. False otherwise.
    '''
    q_score_sum = 0
    for base in seq_quality:
        q_score_sum += ord(base) - 33
    mean_q_score = q_score_sum / len(seq_quality)
    return mean_q_score >= quality_threshold
