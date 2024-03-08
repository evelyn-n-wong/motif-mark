## OOP Motif Mark ##

Alternative splicing is an essential mechanism for producing a diverse array of proteins for gene expression, cellular differentiation, and development. Splicing refers to the process of taking pre-mRNA and removing introns and retaining exons. Multiple transcripts can be made from a single gene locus with different exon-intron junctions pieced together which leads to functions proteins or non-functional mRNA.

- The goal of this program is a python script that is capable of plotting binding motifs on the sequence. 

## Program Description ##

This program contains a python script and a bash script that the user can use to run python script easily. The python script asks the user to input two files: a fasta file and a file of motifs. 

To run the python script directly on your terminal, enter: python /path/to/file/motif_mark-oop.py -f <fasta file name> -m <motifs file name>

The user can also just run the bash script and make adjustments to where the files are located on their local computer. 

- The user should have create a new environment on their local computer and install cairo. Please visit the cairo website for greater details on how to install the cairo packages onto their environment. 

The program outputs a png using the same naming convention as the fasta file that was passed. To see an example of what the image looks like, see Figure_1.png. 

The header contains the name of the FASTA file followed by all identified gene names for each iteration. For each gene, there is a singular line representing the gene length and motifs identified and colored accordingly. Exon on each gene is also identified (as per case conventions (upper_case) in the file). Upstream introns and downstream introns are also colored accordingly. A legend is provided at the bottom to indicate all identified motifs. 

Please see below 'Minimum Requirements' for the requirements assigned to this project.

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