from typing import Union, Tuple

def filter_gc(seq:str, gc_bounds:Union[float, Tuple[float]]) -> bool:
    if isinstance(gc_bounds, Union[float, int]):
        gc_bounds = (0, gc_bounds)
    print(gc_bounds)
    seq = seq.upper()
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    return gc_bounds[0] <= gc_content * 100 <= gc_bounds[1]


print(filter_gc('GGGGG', 80))