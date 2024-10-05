dna_bases = {'A', 'T', 'G', 'C'}
rna_bases = {'A', 'U', 'G', 'C'}
nucleotides = {'A', 'T', 'U', 'G', 'C'}
dna_base_pairs = {'a' : 't', 't' : 'a', 'g' : 'c', 'c' : 'g', 'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
rna_base_pairs = {'a' : 'u', 'u' : 'a', 'g' : 'c', 'c' : 'g', 'A' : 'U', 'U' : 'A', 'G' : 'C', 'C' : 'G'}
codons = {'UUU': 'F', 'UUC': 'F', 
          'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
          'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 
          'UAU': 'Y', 'UAC': 'Y', 
          'UAA': '*', 'UAG': '*', 'UGA': '*', 
          'UGU': 'C', 'UGC': 'C',
          'UGG': 'W', 
          'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
          'CAU': 'H', 'CAC': 'H', 
          'CAA': 'Q', 'CAG': 'Q', 
          'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 
          'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 
          'AUG': 'M!', 
          'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
          'AAU': 'N', 'AAC': 'N', 
          'AAA': 'K', 'AAG': 'K', 
          'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 
          'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
          'GAU': 'D', 'GAC': 'D', 
          'GAA': 'E', 'GAG': 'E', 
          'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
