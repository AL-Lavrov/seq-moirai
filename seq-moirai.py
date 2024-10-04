import dna_rna_tools as drt

def run_dna_rna_tools(*seqs_and_action, **options):
    seqs, action, options = drt.orginize_data(*seqs_and_action, **options)
    print(options)
    return [drt.action_map[action](seq, **options) for seq in seqs]


print(run_dna_rna_tools('GAAAAGGGgc', 'uuuaTaaa', 'complement', fragme=1))