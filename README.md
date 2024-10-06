# Seq-moirai

<img align="right" src="seq_moirai.png" alt="seq-moirai" width="310">

**Seq-moirai** is Python package for simple processing of nucleic acids sequences and basic filtering of FASTQ sequence files.

**Programming language:** Python3

**OS:** Windows, Linux, MacOS

**External dependencies:** None

**License:** CC0

**Version:** 0.2 (06.10.2024)

**Author:** Aleksandr Lavrov


## Contents

- [Seq-moirai](#seq-moirai)
  - [Contents](#contents)
  - [Installation](#installation)
  - [Modules](#modules)
    - [dna\_rna\_tools](#dna_rna_tools)
    - [filter\_fastq](#filter_fastq)
  - [Examples of use](#examples-of-use)
  - [Versions](#versions)
    - [v0.2 (06.10.2024)](#v02-06102024)
    - [v0.1 (29.09.2024)](#v01-29092024)

## Installation

One can get seq-moirai directly from repository:

```
git clone git@github.com:AL-Lavrov/seq-moirai.git
```

Import whole package or particular module needed into your script:

```python
import seq_moirai
# OR
from seq_moirai import run_dna_rna_tools
```

## Modules

seq-moirai contains ~~three~~ two distinct modules.

### dna_rna_tools

Module can be called by `run_dna_rna_tools()` function. This function accepts any number of nucleic acids sequences followed by the desired transformation and optional keyword arguments. No mixed input of DNA/RNA allowed.

**Available actions:**

- `'reverse'` - reverses sequence in 3' to 5' direction;
- `'transcribe'` - transforms DNA sequence;
- `'complement'` - returns complementary sequence in 3'-5' direction;
- `'reverse-complement'` - returns complementary sequence in 5'-3' direction;
- `'count_gc'` - returns a string with GC fraction in sequence;
- `'translate'` - translates the entire sequence into a protein. Start codons are marked as `M!`, and stop codons as `*`. The reading frame is determined by the `frame` argument. Bases not forming a complete triplet are truncated.

**Optional keyword arguments:**

- `seq_type` - must be given as `frame =  'rna'` in cases when RNA input sequences contain no U bases. Overrides mixed RNA/DNA sequence input control.
- `frame` - can take values `1`, `2`, `3`. By default is considered to be `1`.

### filter_fastq

`filter_fastq()` accepts fastq file converted to dictionary (see [examples](#examples-of-use)) and returns a dictionary of the same structure without entries that do not match bounds. Three additional arguments can be provided:

- `gc_bounds` - defines the acceptable GC-fraction range. Sequences outside this range are discarded. Default value is  `(0, 80)`, If a single number is given, it is taken as an upper bound.

- `length_bounds` - defines the acceptable range for sequence length. Sequences that are not included in the interval are discarded. Default value is  `(0, 2**32)`, If a single number is given, it is taken as an upper bound.

- `quality_threshols` - the lower limit of the mean quality score (Phred33). Sequences below this threshold are discarded. Default is `0`.

## Examples of use

1. Translation of protein
   
   ```python
   result = run_dna_rna_tools('uuaugaaucacc', 'gauuaugc', 'auuagccaccgag', 
                              'translate', frame=2)
   print(result)
   # OUTPUT
   # ['YES', 'IM!', 'LATE']
   ```

2. Filtering fastq
   
   See `example_input`, `example_output` in `example` folder.

## Versions

### v0.2 (06.10.2024)

- Added functions translate() and count_gc() are added to dna_rna_tools
- All functions in dna_rna_tools now use same structure of input
- Added filter_fastq module and corresponding function into seq-moirai
- Updated function typing and added docstrings
- Added README.md


### v0.1 (29.09.2024)

- Initial version

