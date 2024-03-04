## OOP Motif Mark ##

## Program Description ##

## Minimum Requirements ##
- Well commented, Python3 compatible, object-oriented code, with CLEAR readme.md file
- Public GitHub repo named motif-mark [WARNING! Assignment will not be graded if repo incorrectly named!]
- Script named motif-mark-oop.py [WARNING! Assignment will not be graded if script incorrectly named!]
- Use argparse with options: [WARNING! Assignment will not be graded if argparse options not as requested!]
    - f: fasta file
    - m: motifs file
- Output file has same prefix as input file (e.g. Figure_1.fa -> Figure_1.png)
- Input FASTA file (seqs ≤1000 bases) and motifs file (≤10 bases each, one motif per line in a text file)
- Capable of handling:
    - Motifs with ambiguous nucleotides (see https://en.wikipedia.org/wiki/Nucleic_acid_notation Links to an external site.)
    - Multiple sequences (max 10 in the data you will be provided)
    - Multiple motifs (max 5 in the data you will be provided)
- Consider:
    - How you will handle overlapping motifs
    - How you will denote introns/exons
- All features (motifs, introns, exons) should be to scale
- Output single, well-labeled figure, per FASTA file
- png output
- Key/labeling
- Question: How are exons and introns defined in the FASTA file? Don't start planning or coding until you're sure you know the correct answer (it's not a trick question).
- Should be able to be run in the following environment:
    - conda create -n my_pycairo pycairo
    - conda activate my_pycairo
    - [WARNING! Assignment will not be graded if other packages are required!]
- Style: In Python, the convention is to use CamelCase for your class names, and underscore / "snake_case" style for your objects. For example your class might be "BioformaticsLover" and an object derived from that class might be "my_happy_bioinformatics_student".  Shame on Jason for not sharing that during the OOP workshop. See this Links to an external site.. We encourage you to do the same.-