import modules.dna_rna_tools as drt
from typing import Union

def run_dna_rna_tools(*seqs_and_action: str, **options: drt.Options) -> Union[str, list]:
    seqs, action, options = drt.get_data(*seqs_and_action, **options)
    action_map = {'transcribe': drt.transcribe,
              'reverse': drt.reverse,
              'reverse_complement': drt.reverse_complement,
              'complement': drt.complement}
    processed_seq = [action_map[action](seq, **options) for seq in seqs]
    return processed_seq[0] if len(processed_seq) == 1 else processed_seq


print(run_dna_rna_tools('aaaccc', 'AAAAAAAA', 'cccccG', 'complement', seq_type=1))