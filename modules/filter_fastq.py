from typing import Union, Tuple


def filter_gc(seq: str, gc_bounds: Union[float, Tuple[float]]) -> bool:
    if isinstance(gc_bounds, Union[float, int]):
        gc_bounds = (0, gc_bounds)
    seq = seq.upper()
    gc_content = (seq.count("G") + seq.count("C")) / len(seq)
    return gc_bounds[0] <= gc_content * 100 <= gc_bounds[1]


def filter_length(seq: str, length_bounds: Union[int, Tuple[int]]) -> bool:
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def filter_quality(seq_quality: str, quality_threshold: int) -> bool:
    q_score_sum = 0
    for base in seq_quality:
        q_score_sum += ord(base) - 33
    mean_q_score = q_score_sum / len(seq_quality)
    return mean_q_score >= quality_threshold
