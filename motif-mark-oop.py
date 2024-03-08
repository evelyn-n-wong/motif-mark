#!/usr/bin/env python
# Author: Evelyn Wong
# Date Created: 2024-02-20
# Last updated: 2024-03-07
# Description: This program script takes a FASTA file and motifs file to plot protein binding motifs on a png image of an exon and flanking introns.
# Collaborators: Tam Ho and Sydney Hamilton

import cairo # cairo package for drawing 
import argparse # get user arguments 
import sys 
import re #regex
import random

class Gene:

    def __init__(self, header: str, seq: str):
        '''Initializing gene by getting FASTA header, sequence, and sequence length'''
        self.header = header # FASTA header
        self.seq = seq # gene sequence 
        self.length = len(seq) # total length of sequence
        self.name = header.split(" ")[0][1:] # gene name from header
        self.chr = header.split(" ")[1] # gene location from header line
        self.exon_seq: list = [] # list to store exon sequences
        self.pos = re.search('\d+-\d+', header).group(0) # gene position on chromosome

    # function to get positions 

    def find_exon_intron(self, seq):
        '''Parse fasta sequence to identify exon start/end positions. Also identifies introns both upstream and downstream'''
        
        # list of tuple (start and ending pos of match of exon in fasta seq)
        exon = [(match.start(), match.end(), match.group()) for match in re.finditer(pattern="[A-Z]+", string=self.seq)]

        introns_upstream = [] # list to hold tuple of starting and ending pos of introns upstream of each exon
        introns_downstream = [] # list to hold tuple of starting and ending pos of introns downstream of each exon
    
        for i in range(len(exon) - 1):
            intron_upstream_end = exon[i][0]  # end pos of upstream intron is start pos of current exon
            intron_downstream_start = exon[i][1]  # start pos of downstream intron is end pos of current exon
            introns_upstream.append((0, intron_upstream_end))  # tuple of intron upstream/downstream pos to list
            introns_downstream.append((intron_downstream_start, exon[i + 1][0]))  # tuple of start pos of downstream intron and start pos of next exon

        # first and last introns
        if exon:
            introns_upstream.insert(0, (0, exon[0][0]))  # start from pos 0 and end at start of first exon[0][0]
            introns_downstream.append((exon[-1][1], len(self.seq)))  # start from end pos of last exon until rest of sequence
            self.exon_seq.extend([exon_seq[2] for exon_seq in exon]) # list of exon sequences without pos

        return exon, introns_upstream, introns_downstream

class Motif:

    # way to handle ambiguous, degenerate ones
    iupac_codes = {
        "A": "[Aa]",
        "T": "[TtUu]",
        "U": "[Tt]",
        "C": "[Cc]",
        "G": "[Gg]",
        "R": "[AaGg]",
        "Y": "[CcTtUu]",
        "W": "[AaTtUu]",
        "S": "[GgCc]",
        "K": "[GgTtUu]",
        "M": "[AaCc]",
        "B": "[CcGgTtUu]",
        "D": "[AaGgTtUu]",
        "H": "[AaCcTtUu]",
        "V": "[AaCcGg]",
        "N": "[AaCcGgTtUu]",
        "Z": "[]"
    }

    def __init__(self, motif_seq, og_motif):
        '''Initialize motif instance'''
        self.motif_seq = motif_seq.upper()
        self.og_motif = og_motif.upper()

    def search_motif(self, sequence):
        '''This function takes in a gene sequence and gets the positions of motifs in the sequence'''
        #print(f'Motif: {self.motif_seq}')
        motif_matches_dict = {}

        for match in re.finditer(self.motif_seq, sequence.upper()):
            start_pos, end_pos = match.start(), match.end()
            motif_matches_dict[(start_pos, end_pos)] = {
                'matched_motif': sequence[start_pos:end_pos],
                'original_motif': self.og_motif
            }
        
        return motif_matches_dict

class BobRoss:
    # class for painting using pycairo

    def __init__(self, file_prefix, gene_obj, max_length, all_gene_motif_pos, exon_gene_dict):
        '''Initalize variables for drawing'''
        self.file_prefix = file_prefix # FASTA file name prefix
        self.max_length = max_length # get max for setting image width
        self.all_gene_motif_pos = all_gene_motif_pos # dictionary of Gene Names and related Motifs + positions
        self.exon_gene_dict = exon_gene_dict # Dictionary containing Key: Gene Name and Value: Tuple of start/end pos for exon
        self.image_height = 800

    def draw_figure_title(ctx, file_prefix):
        ctx.move_to(30,30)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(24)
        font_face = cairo.ToyFontFace("sans-serif", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        ctx.set_font_face(font_face)
        ctx.show_text(f'FASTA file: {file_prefix}')
        
    def draw_gene(ctx, gene_name, gene_len, line_y):

        # draw horizonatal line representing gene sequence length
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(5)
        ctx.move_to(50, line_y)
        ctx.line_to(50 + gene_len, line_y)
        ctx.stroke()
        
        # label for gene names
        ctx.move_to(50, line_y - 35)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(18)
        ctx.show_text(gene_name)
        #print(f'{all_gene_motif_pos[0]} {gene_obj.chrom} {gene_obj.pos}')

    def get_palette(motif_list):
        '''This function returns a list of palette colors depending on the length of motif_list passed in'''
        motif_count = 0 # counter for motifs

        full_palette = [[0, 0, 1],
                        [1, 0, 0],
                        [0, 1, 0],
                        [0.4, 0.4, 0.4],
                        [0.3, 0.3, 0.3]]

        choose_colors = full_palette[0:len(motif_list)]

        return(choose_colors)

    def get_motif_colors(dictionary, color_list):
        '''This function returns dictionary where key = motif, and value = color from the passed_in list of colors'''
        motif_color = {} # holds key of motifs and values of color
        i = 0 # counter for the color_palette

        for motif_info in dictionary.values():
            for motif, pos in motif_info.items():
                if motif not in motif_color:
                    motif_color[motif] = color_list[i]
                    i+=1
        
        return(motif_color) # e.g. {'YGCY': [0.5, 0.5, 0.5], 'GCAUG': [0.2, 0.2, 0.2], 'CATAG': [0.7, 0.7, 0.7], 'YYYYYYYYYY': [0.4, 0.4, 0.4]}
    
    def draw_motif(ctx, motif, start, end, color_assign:dict, line_y):
        '''This function takes in the start, end positions from the all_gene_motif_pos dictionary
        in the main drawing gene function.
        Returns: drawing for motif'''
    
        temp_color = color_assign.get(motif)
        r, g, b = temp_color[0], temp_color[1], temp_color[2] # gets the assigned color for specific motif
        ctx.set_source_rgb(r, g, b)
        ctx.move_to(20 + start, line_y)
        ctx.set_line_width(35)
        ctx.line_to(20 + end, line_y)
        ctx.stroke()


    def draw_exon(ctx, start, end, line_y):
        '''This function takes in the context for drawing and starting, ending positions from the exon_gene_dict.
        exon_gene_dict: key is the gene_name, value is a tuple of the start/end position of the xon on that gene 
            e.g. Gene: INSR, Pos: (276, 312) 
        Returns: drawing for exon per gene'''
        print(f'Exon_pos: Start:{start}, End: {end}')
        ctx.set_source_rgba(1, 0, 1, 0.1) # set color
        ctx.set_line_width(50)
        ctx.move_to(20 + start, line_y) # fix
        ctx.line_to(20 + end, line_y) # fix
        ctx.stroke()
        print("Exon drawn successfully")

    def draw_ds_intron(ctx, start, end, line_y):
        '''This function takes in the context for drawing and starting, ending positions from the exon_gene_dict.
        exon_gene_dict: key is the gene_name, value is a tuple of the start/end position of the xon on that gene 
            e.g. Gene: INSR, Pos: (276, 312) 
        Returns: drawing for exon per gene'''
        print(f'Exon_pos: Start:{start}, End: {end}')
        ctx.set_source_rgb(1, 1, 0) # set color
        ctx.set_line_width(10)
        ctx.move_to(20 + start, line_y) # fix
        ctx.line_to(20 + end, line_y) # fix
        ctx.stroke()
        print("Downstream intron drawn successfully")
    
    def draw_up_intron(ctx, start, end, line_y):
        '''This function takes in the context for drawing and starting, ending positions from the exon_gene_dict.
        exon_gene_dict: key is the gene_name, value is a tuple of the start/end position of the xon on that gene 
            e.g. Gene: INSR, Pos: (276, 312) 
        Returns: drawing for exon per gene'''
        print(f'Exon_pos: Start:{start}, End: {end}')
        ctx.set_source_rgb(0, 1, 1) # set color
        ctx.set_line_width(10)
        ctx.move_to(20 + start, line_y) # fix
        ctx.line_to(20 + end, line_y) # fix
        ctx.stroke()
        print("Upstream intron drawn successfully")

    def draw_legend(ctx, color_assign: dict, legend_x, legend_y, legend_font_size, legend_spacing):
        ''' This function takes in coordinates for drawing a legend ''' 
        for motif, color in color_assign.items():
            ctx.set_source_rgb(*color)
            ctx.rectangle(legend_x, legend_y, 10, 10)
            ctx.fill()
            ctx.move_to(legend_x + 20, legend_y + 10)
            ctx.set_font_size(legend_font_size)
            ctx.show_text(motif)
            legend_y += legend_spacing 

def oneline_fasta(fasta_file):
    '''This function takes in a multiline FASTA file and returns a one-line fasta'''

    fasta_dict = {} # dictionary to hold fasta header as key and one-line sequence as value
    current_header = ""

    with open(fasta_file, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                current_header = line.strip("\n")
                fasta_dict[current_header] = "" # add header  
            else:
                fasta_dict[current_header] += line.strip("\n")

    return fasta_dict

def convert_motif(motif_file):
    '''This function takes in a motif file and returns a list of motifs, all upper-case converted using iupac codes for generate cases.
    Returns list of converted str: motif sequences'''
    motif_list: list = [] # holds regex for converted motifs
    original_list: list = [] # holds original motifs 
    with open(motif_file, "r") as mf:
        for line in mf:
            new_motif = ""
            line = line.upper().strip() # strip line; original motif
            for b in line:
                new_motif += Motif.iupac_codes.get(b,b)
            motif_list.append(new_motif)
            original_list.append(line)
    return motif_list, original_list

if __name__ == '__main__':
    ## Main Function

    def get_args():
        '''Get files from user. See 'motif_mark.sh' to enter file information easily.'''
        parser = argparse.ArgumentParser(description="This function takes in user arguments for demultiplexing. Arguments required: -f, -m, -o")
        parser.add_argument("-f", "--file", help = "designate absolute file path to fasta file", required = True)
        parser.add_argument("-m", "--motif", help = "designates absolute file path to motifs file", required = True) #e.g. "Figure_1a.fa" -> "Figure_1.png" 
        args = parser.parse_args()
        return args

    # Setting up variables 
    args = get_args()
    fasta_file = args.file
    motif_file = args.motif

    # get converted motif_list
    motif_list, original_list = convert_motif(motif_file) # list of converted regex motifs and original motifs
    # get one line FASTA (key: headers, values: sequences)
    fasta_dict = oneline_fasta(fasta_file)
    # extract file prefix from fasta file for image output
    split_fasta = re.split("\.", fasta_file)
    file_prefix = split_fasta[0] # use for output
    output_file_path = f'{file_prefix}.png'
    # instantiate motif objects 
    motif_objects = [Motif(motif_seq, original_motif) for motif_seq, original_motif in zip(motif_list, original_list)]
    # instantiate gene objects
    gene_objs = [Gene(header, sequence) for header, sequence in fasta_dict.items()]
    # get sequences for all genes in a list
    gene_seqs = [values for key, values in fasta_dict.items()]

    all_gene_motif_pos = {} # dictionary to hold all matched positions for each gene
    exon_gene_dict = {} # hold gene name as key and exon positions as values
    intron_up_dict = {} # hold gene name as key and intron_upstream positions as values
    intron_down_dict = {} # hold gene name as key and intron_downstream_positions as values 
    max_length: int = 0 # length of max_sequence in gene_seqs
    sequence_length_list = [len(gene_seq) for gene_seq in gene_seqs] # list of gene sequence lengths 

    # Process each gene
    for gene_obj in gene_objs:
        # get max length of sequence 
        if gene_obj.length > max_length:
            max_length = gene_obj.length
        # Find motifs in the gene sequence 
        gene_motif_positions = {} # dictionary to hold all matched positions for current gene object
        # print(gene_obj.name)
        for motif_obj in motif_objects:
            motif_matches_dict = motif_obj.search_motif(gene_obj.seq) # returns dictionary; key: tuple of position, value: tuple of matched motif and og motif
            #print(f'Motif Matches Dict Pos: {motif_matches_dict}')

            # add matched positions of each motif to new dictionary 
            for pos_info, motif_info in motif_matches_dict.items():
                original_motif = motif_info['original_motif']
                if original_motif not in gene_motif_positions:
                    gene_motif_positions[original_motif] = []
                gene_motif_positions[original_motif].append(pos_info)

        # instantiate exons, upstream introns, and downstream introns
        exon, intron_up, intron_down = gene_obj.find_exon_intron(gene_obj.seq)
        # print(f'Gene Length: {gene_obj.length}')
        # print(f'Exon: {exon}')
        # extend found positions with exon 
        for exon_start, exon_end, exon_seq in exon:
            for motif_name, motif_positions in gene_motif_positions.items():
                for i, motif_pos in enumerate(motif_positions):
                    start, end = motif_pos
                    if exon_start <= start <= exon_end:
                        # update start pos based on exon positions to get all motif positions for each gene
                        gene_motif_positions[motif_name][i] = (start + exon_start, end + exon_end)
            exon_gene_dict[gene_obj.name] = (exon_start, exon_end)
        
        for intron_up_start, intron_up_end in intron_up:
            intron_up_dict[gene_obj.name] = (intron_up_start, intron_up_end)

        for intron_down_start, intron_down_end in intron_down:
            intron_down_dict[gene_obj.name] = (intron_down_start, intron_down_end)
                        

        # add gene_motif_pos to all_gene_motif_pos dictionary
        all_gene_motif_pos[gene_obj.name] = gene_motif_positions

    # getting only singular motifs as a list
    all_motif_list = []
    for motif_dict in all_gene_motif_pos.values():
        for motifs in motif_dict.keys():
            if motifs not in all_motif_list:
                all_motif_list.append(motifs)

    for gene, intron_up in intron_up_dict.items():
        print(f'Gene: {gene}, Pos: {intron_up}\n')

    for gene, intron_down in intron_down_dict.items():
        print(f'Gene: {gene}, Pos: {intron_down}\n')

    for gene, exon in exon_gene_dict.items():
        print(f'Gene: {gene}, Pos: {exon}\n')
    # for gene_obj in gene_objs:
    #     #Draw the gene
    #     print(gene_obj.name)
    #     br = BobRoss(file_prefix, gene_obj, max_length, all_gene_motif_pos, exon_gene_dict)
    

    ## Setting up Drawing function
    image_width = 1000
    image_height = 600

    sfc = cairo.ImageSurface(cairo.FORMAT_ARGB32, image_width, image_height)
    ctx = cairo.Context(sfc)
    ctx.set_source_rgb(1,1,1)
    ctx.paint()
    
    c_list = BobRoss.get_palette(all_motif_list) # gets a list of color based on number of motifs
    
    color_assign = BobRoss.get_motif_colors(all_gene_motif_pos, c_list) # assigns colors to motifs found across genes
    #e.g. ATP2A1, {'YGCY': [(31, 35), (40, 44), (45, 49), (97, 101), (113, 117), (148, 152), (183, 187), (194, 198), (208, 212), (523, 654), (538, 669), (543, 674), (549, 680), (558, 689), (366, 370), (391, 395), (397, 401), (403, 407), (408, 412), (428, 432), (475, 479), (495, 499), (500, 504), (510, 514), (521, 525), (535, 539)], 'CATAG': [(228, 233)], 'YYYYYYYYYY': [(17, 27), (474, 611), (514, 651)]}

    s_counter = 0 

    l = 100
    for gene, intron_down_pos in intron_down_dict.items():
        print(f'Gene {gene} Downstream intron pos {intron_down_pos}')
        i_down_start, i_down_end = intron_down_pos[0], intron_down_pos[1]
        BobRoss.draw_ds_intron(ctx = ctx, start = i_down_start, end = i_down_end, line_y = l)
        l += 100
    
    k = 100

    for gene, intron_up_pos in intron_up_dict.items():
        print(f'Gene {gene} Upstream intron pos: {intron_up_pos}')
        i_up_start, i_up_end = intron_up_pos[0], intron_up_pos[1]
        BobRoss.draw_up_intron(ctx = ctx, start = i_up_start, end = i_up_end, line_y = k)
        gene_l = sequence_length_list[s_counter]
        BobRoss.draw_gene(ctx = ctx, gene_name = gene, gene_len = gene_l, line_y = k)
        k += 100

    j = 100
    for gene, exon_pos in exon_gene_dict.items():
        #print(f'Gene {gene} Exon_pos, {exon_pos}')
        exon_start, exon_end = exon_pos[0], exon_pos[1]
        BobRoss.draw_exon(ctx = ctx, start = exon_start, end = exon_end, line_y = j)
        j += 100

    i = 100
    for motif_info in all_gene_motif_pos.values():
        #print("Motif_info", motif_info)
        for motif, pos in motif_info.items():
            #print("Full:", pos)
            for x, y in pos:
                current_start, current_end = x, y
                BobRoss.draw_motif(ctx = ctx, motif = motif, start = current_start, end = current_end, color_assign = color_assign, line_y = i)
        i += 100

    BobRoss.draw_figure_title(ctx, file_prefix)

    for c in c_list:
        print(c)

    legend_x = 50
    legend_y = 450
    legend_spacing = 20
    legend_font_size = 12

    BobRoss.draw_legend( ctx = ctx, color_assign = color_assign, legend_x = legend_x, legend_y = legend_y, legend_font_size = legend_font_size, legend_spacing = legend_spacing)

    ## Output saved image as png 

    sfc.write_to_png(output_file_path)
    print(f'Image saved to: {output_file_path}')