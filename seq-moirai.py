import modules.dna_rna_tools as drt
from typing import Union, Tuple, Dict


def run_dna_rna_tools(*seqs_and_action: str, **options: drt.Options
                      ) -> Union[str, list]:
    
    seqs, action, options = drt.get_data(*seqs_and_action, **options)
    action_map = {'transcribe': drt.transcribe,
                  'reverse': drt.reverse,
                  'reverse_complement': drt.reverse_complement,
                  'complement': drt.complement,
                  'count_gc': drt.count_gc,
                  'translate': drt.translate}
    processed_seq = [action_map[action](seq, **options) for seq in seqs]
    return processed_seq[0] if len(processed_seq) == 1 else processed_seq


def filter_fastq(seqs: Dict[str, Tuple[str]], 
                 gc_bounds: Union[float, Tuple[float]] = (0, 100), 
                 length_bounds: Union[float, Tuple[float]] = (0, 2**32), 
                 quality_bounds: float = 0
                 ) -> Dict[str, Tuple[str]]:
    pass