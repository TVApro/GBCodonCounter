#! /usr/bin/python3
# Программа для получения статистических данных о кодонах в геноме

'''The instruction from the Codon_Counter.py file

Analysis requires pandas, openpyxl and collections modules:

 pip3 install pandas
 pip3 install openpyxl
 pip3 install collections

The program accepts as input files:

*.gb, .gbk, .gbff and other GenBank files containing annotations

*.fna files containing nucleotide sequences of genes (be careful - the program is not yet able to adequately exclude tRNA, rRNA and ncRNA from such files)

The program generates *.fna files without RNA itself in a separate folder, they can be used later
P.S.

Warning 1. The files will contain numbers with dots, which are not defined as numbers in Russian localization.

Solution 1: "Edit - Replace All" in any office program (MS Excel, OnlyOffice, LibreOffice, WPS Office etc.), replace all dots with commas

Warning 2. The error of determined percentages is set to the 6th decimal place by default. Since the average number of codons in a bacterial genome is approximately 600,000 - 700,000, for eukaryotic comparisons or viruses someone may need to expand or narrow the margin of error.

Solution 2: Line 279: X = format(float((X/summ)*100), '.6f') Instead of 6, put any desired number of decimal places

Happy exploration!'''

import re
from re import sub
import textwrap as tw
import pandas as pd
import openpyxl
import os
from collections import Counter

def start():
    x1 = False
    global single_mode
    global group_mode
    global percentage
    print('\nWelcome to the codon counting program!')
    while x1 == False:
        q = input('Select operating mode: \n 1 - single analysis \n 2 - group analysis \n 3 - exit \n\n Answer: ')
        otvet_list = ('1', '2', '3')
        if q == '1': 
            single_mode = True
            group_mode = False
            print('Single genome analysis mode\n')
            x1 = True
        if q == '2': 
            group_mode = True
            single_mode = False
            print('Group of genomes analysis mode\n')
            x1 = True
        if q == '3':
            print('Exit')
            exit()
        if q not in otvet_list: 
            print('Invalid option selected. Start selection again\n')
    w = input('\nNeed a percentage?\n 0 - no\n pass - default yes\n\n Answer: ')
    if w == '0':
        percentage = False
        print('Determination of the absolute number of amino acids is chosen\n')
    else:
        print('Percentage calculation selected\n')


def open_file(seq):
    with open(seq, 'r') as file:
            text = file.read()
    print('Opened', file)
    return text
   
def file_format_def(text):
    forma = ''
    if text[:1] == '>':
        forma = 'fasta'
        print('File formate: fasta')
    if text[:5] == "LOCUS":
        forma = 'gbk'
        print('File formate: genbank')
    if text[:1] != '>' and text[:5] != "LOCUS":
        forma = 'error'
        print('The input file does not match the format (genbank or .fna)')
    return forma
        
def gbk_to_fna(text, file_name, path):
    # Регулярные выражения
    global all_AA_count
    gene_name_reg = re.compile(r'/product=\".*?\"', re.DOTALL)
    gene_id_reg = re.compile(r'/protein_id=\".*?\"', re.DOTALL)
    annot_split_reg = re.compile(r'/.+=\\')
    LOCUS = re.compile(r'LOCUS   .*')
    title = re.compile(r' {5}\S+  +')
    edges = re.compile(r'\d+\W+\d+|complement\(\d+\W+\d+?\)')
    
    x = 0
    prot_number = 0
    text_list = []

    def beautifull_sequence(seq):
        seq = tw.fill(seq, width=60)
        return seq

    def simple_sequence(seq):
        seq = seq.replace(' ', '').replace('\n', '').replace('-', '')
        seq = seq.replace('a', 'A').replace('c', 'C').replace('t', 'T').replace('g', 'G').replace('n', 'N')
        return seq

    def complement_sequence(seq):
        seq = seq.replace(' ', '').replace('\n', '').replace('-', '')
        seq = seq.replace('A', '%1%').replace('a', '%1%').replace('T', 'A').replace('t', 'A').replace('%1%', 'T')
        seq = seq.replace('C', '%2%').replace('c', '%2%').replace('G', 'C').replace('g', 'C').replace('%2%', 'G')
        seq = seq[::-1]
        return seq

    text_list = re.split(LOCUS, text) # разбиваем геном на участки LOCUS
    for i in text_list:
        if i:
            # поиск нуклеотидной последовательности
            ind_FEATURES = i.find("FEATURES ")
            ind_ORIGIN = i.find("ORIGIN")
            ind_end = i.find("//", ind_ORIGIN)
            seq =  simple_sequence(re.sub(r'\d*', '', i[ind_ORIGIN+6:ind_end]))
            annotations = i[ind_FEATURES:ind_ORIGIN]
            gene_list = re.split(title, annotations)
            for j in gene_list:
                    if '/product="' in j:
                        if '/translation=' in j:
                            gene_name_position = re.search(gene_name_reg, j)
                            gene_name = re.sub(r' +', ' ', (gene_name_position.group()).replace('\n', ''))
                            g_edges = re.search(edges, j)
                            gg_edges = g_edges.group()
                            ggg_edges = re.findall('\d+', gg_edges)
                            u = 0
                            na_sequence = str('')
                            for h in ggg_edges:
                                u += 1
                                if u == 1:
                                    st = int(h)
                                if u == 2:
                                    fin = int(h)
                                    delta_l = abs(fin-st)
                                    na_sequence = seq[st-1:st+delta_l]
                                    break
                            if 'complement' in gg_edges:
                                na_sequence = complement_sequence(na_sequence)
                                        
                            prot_number += 1

                            fna_name = gene_name.replace('/product=\"', '>').replace('"', '')
                            nucleic_acids = beautifull_sequence(na_sequence)
                            if gene_name == '>' or gene_name == '> ':
                                gene_name = '>Unknown_protein_' + str(lambda x: x+1)
                            
                            fnafilename = file_name+'.fna'
                            new_path = os.path.join(path, fnafilename)
                            if prot_number == 1:
                                with open(new_path, 'w') as fna_file:
                                    print(fna_name+'\n'+nucleic_acids, file=fna_file)
                            else:
                                with open(new_path, 'a') as fna_file:
                                    print(fna_name+'\n'+nucleic_acids, file=fna_file)
    all_AA_count += prot_number
    print('Read', prot_number, 'protein sequences')

def fna_in_nucleic_counter(file):
    global all_codon_count
    def codon_obtain(seq):
        ind_codon_list = []
        krat = len(seq) // 3 # целочисленное деление
        times = 1
        while times <= krat:
            for x in range(0, len(seq), 3):
                codon = seq[x : 3 + x] # не учитываем последние нуклеотиды, если они меньше кодона
                ind_codon_list.append(codon)
                times += 1
        return ind_codon_list

    A_new = 0
    G_new = 0
    T_new = 0
    C_new = 0

    with open(file, 'r') as myfile:
        text = myfile.read()
    ALL_CODONS = []
    fasta_list = text.split(">")

    for i in fasta_list:
        if i:
            end = i.find("\n")
            fasta_name = '>'+sub('\n', '', i[:end])
            fasta_sequence = sub('\n', '', i[end:])
            A_new += fasta_sequence.count('A')
            G_new += fasta_sequence.count('G')
            T_new += fasta_sequence.count('T')
            C_new += fasta_sequence.count('C')
            fasta_codons = codon_obtain(fasta_sequence)
            ALL_CODONS = ALL_CODONS + fasta_codons # полный список всех кодонов генома
    genome_name_ind = file.rfind('/')
    genome_name = file[genome_name_ind+1:]
    all_codon_count += len(ALL_CODONS)
    GC = (G_new+C_new)*100/(G_new+C_new+A_new+T_new)
    perGC = str(format(GC, '.2f'))+'%'
    print('GC-composition of the coding part', perGC)

    CODONS_list = Counter(ALL_CODONS)
    
    TGT = CODONS_list['TGT'] 
    TGC = CODONS_list['TGC']
    TGG = CODONS_list['TGG']
    GAT = CODONS_list['GAT']
    GAC = CODONS_list['GAC']
    TTT = CODONS_list['TTT']
    TTC = CODONS_list['TTC']
    GGT = CODONS_list['GGT']
    GGC = CODONS_list['GGC']
    GGA = CODONS_list['GGA']
    GGG = CODONS_list['GGG']
    ACT = CODONS_list['ACT']
    ACC = CODONS_list['ACC']
    ACA = CODONS_list['ACA']
    ACG = CODONS_list['ACG']
    TCT = CODONS_list['TCT']
    TCC = CODONS_list['TCC']
    TCA = CODONS_list['TCA']
    TCG = CODONS_list['TCG']
    AGT = CODONS_list['AGT']
    AGC = CODONS_list['AGC']
    ATG = CODONS_list['ATG']
    GCT = CODONS_list['GCT']
    GCC = CODONS_list['GCC']
    GCA = CODONS_list['GCA']
    GCG = CODONS_list['GCG']
    TAT = CODONS_list['TAT']
    TAC = CODONS_list['TAC']
    CAT = CODONS_list['CAT']
    CAC = CODONS_list['CAC']
    CTT = CODONS_list['CTT']
    CTC = CODONS_list['CTC']
    CTA = CODONS_list['CTA']
    CTG = CODONS_list['CTG']
    TTA = CODONS_list['TTA']
    TTG = CODONS_list['TTG']
    GAA = CODONS_list['GAA']
    GAG = CODONS_list['GAG']
    CCT = CODONS_list['CCT']
    CCC = CODONS_list['CCC']
    CCA = CODONS_list['CCA']
    CCG = CODONS_list['CCG']
    GTT = CODONS_list['GTT']
    GTC = CODONS_list['GTC']
    GTA = CODONS_list['GTA']
    GTG = CODONS_list['GTG']
    CGT = CODONS_list['CGT']
    CGC = CODONS_list['CGC']
    CGA = CODONS_list['CGA']
    CGG = CODONS_list['CGG']
    AGA = CODONS_list['AGA']
    AGG = CODONS_list['AGG']
    AAA = CODONS_list['AAA']
    AAG = CODONS_list['AAG']
    AAT = CODONS_list['AAT']
    AAC = CODONS_list['AAC']
    CAA = CODONS_list['CAA']
    CAG = CODONS_list['CAG']
    ATT = CODONS_list['ATT']
    ATC = CODONS_list['ATC']
    ATA = CODONS_list['ATA']
    TAA = CODONS_list['TAA']
    TGA = CODONS_list['TGA']
    TAG = CODONS_list['TAG']

    new_list = [TGT, TGC, TGG, GAT, GAC, TTT, TTC, GGT, GGC, 
                GGA, GGG, ACT, ACC, ACA, ACG, TCT, TCC, TCA, TCG,
                AGT, AGC, ATG, GCT, GCC, GCA, GCG, TAT, TAC, CAT,
                CAC, CTT, CTC, CTA, CTG, TTA, TTG, GAA, GAG, CCT,
                CCC, CCA, CCG, GTT, GTC, GTA, GTG, CGT, CGC, CGA,
                CGG, AGA, AGG, AAA, AAG, AAT, AAC, CAA, CAG, ATT,
                ATC, ATA, TAA, TGA, TAG]

    global percentage
    if percentage == False:
        new_row = pd.DataFrame({genome_name: new_list})

    if percentage == True:
        def percentof(X, summ):
            X = format(float((X/summ)*100), '.6f')
            return X

        new_per_list = []
        for i in new_list:
            i_new = percentof(i, len(ALL_CODONS))
            new_per_list.append(i_new) 

        new_row = pd.DataFrame(new_per_list, columns = [genome_name])
    return new_row

# начало программы
single_mode = False
group_mode = False
percentage = True
new_row = ()
all_codon_count = 0
all_AA_count = 0

first_column = pd.DataFrame({'AA': ['Cys', 'Cys', 
                                    'Trp', 
                                    'Asp', 'Asp', 
                                    'Phe', 'Phe', 
                                    'Gly', 'Gly', 'Gly', 'Gly', 
                                    'Thr', 'Thr', 'Thr', 'Thr',
                                    'Ser', 'Ser', 'Ser', 'Ser', 'Ser', 'Ser', 
                                    'Met', 
                                    'Ala', 'Ala', 'Ala', 'Ala', 
                                    'Tyr', 'Tyr', 
                                    'His', 'His', 
                                    'Leu', 'Leu', 'Leu', 'Leu', 'Leu', 'Leu', 
                                    'Glu', 'Glu', 
                                    'Pro', 'Pro', 'Pro', 'Pro',
                                    'Val', 'Val', 'Val', 'Val', 
                                    'Arg', 'Arg', 'Arg', 'Arg', 'Arg', 'Arg', 
                                    'Lys', 'Lys', 
                                    'Asn', 'Asn', 
                                    'Gln', 'Gln', 
                                    'Ile', 'Ile', 'Ile', 
                                    'STOP', 'STOP', 'STOP']})

Cys_list = pd.DataFrame(['TGT', 'TGC'], columns=['CODONS'])
Trp_list = pd.DataFrame(['TGG'], columns=['CODONS'])
Asp_list = pd.DataFrame(['GAT', 'GAC'], columns=['CODONS'])
Phe_list = pd.DataFrame(['TTT', 'TTC'], columns=['CODONS'])
Gly_list = pd.DataFrame(['GGT', 'GGC', 'GGA', 'GGG'], columns=['CODONS'])
Thr_list = pd.DataFrame(['ACT', 'ACC', 'ACA', 'ACG'], columns=['CODONS'])
Ser_list = pd.DataFrame(['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], columns=['CODONS'])
Met_list = pd.DataFrame(['ATG'], columns=['CODONS'])
Ala_list = pd.DataFrame(['GCT', 'GCC', 'GCA', 'GCG'], columns=['CODONS'])
Tyr_list = pd.DataFrame(['TAT', 'TAC'], columns=['CODONS'])
His_list = pd.DataFrame(['CAT', 'CAC'], columns=['CODONS'])
Leu_list = pd.DataFrame(['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'], columns=['CODONS'])
Glu_list = pd.DataFrame(['GAA', 'GAG'], columns=['CODONS'])
Pro_list = pd.DataFrame(['CCT', 'CCC', 'CCA', 'CCG'], columns=['CODONS'])
Val_list = pd.DataFrame(['GTT', 'GTC', 'GTA', 'GTG'], columns=['CODONS'])
Arg_list = pd.DataFrame(['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], columns=['CODONS'])
Lys_list = pd.DataFrame(['AAA', 'AAG'], columns=['CODONS'])
Asn_list = pd.DataFrame(['AAT', 'AAC'], columns=['CODONS'])
Gln_list = pd.DataFrame(['CAA', 'CAG'], columns=['CODONS'])
Ile_list = pd.DataFrame(['ATT', 'ATC', 'ATA'], columns=['CODONS'])
STOP_list = pd.DataFrame(['TAA', 'TGA', 'TAG'], columns=['CODONS'])

new_frame = pd.concat([Cys_list, Trp_list, Asp_list, Phe_list,
                    Gly_list, Thr_list, Ser_list, Met_list, Ala_list,
                    Tyr_list, His_list, Leu_list, Glu_list, Pro_list,
                    Val_list, Arg_list, Lys_list, Asn_list, Gln_list,
                    Ile_list, STOP_list], axis=0)
new_frame.reset_index(drop=True, inplace=True)
start_frame = pd.concat([first_column, new_frame], axis=1)

start()
if single_mode == True:
    genome_name_ind = file.rfind('/')
    path = file[:genome_name_ind-1]
    try:
        x2 = False
        while x2 == False:
            file = input('Enter sequence address: \n')
            genome=open_file(file)
            forma = file_format_def(genome)
            if forma == 'error':
                print('\nTry again\n')
            else:
                x2 = True
        if forma == 'gbk':
            output_file = input('Select the name of the table file with the analysis results: ')
            print('\n')
            gbk_to_fna(genome, file, path)
            newfile = file+'.fna'
            new_row = fna_in_nucleic_counter(newfile)
        if forma == 'fasta':
            output_file = input('Select the name of the table file with the analysis results: ')
            print('\n')
            new_row = fna_in_nucleic_counter(file)
        start_frame = pd.concat([start_frame, new_row], axis=1)
        start_frame.to_excel(output_file+'.xlsx')
        print('\nRecorded file '+output_file+'.xlsx')
        print('\nTotal processed', all_codon_count, 'codons')
    except:
        print('Failed to write file')


if group_mode == True:
    genomes_list = []
    error_number = 0
    genome_number = 0
    x3 = False
    error_format = False
    while x3 == False:
        path = input("Specify the path to the folder with genomes: ")
        if os.path.isdir(path):
            genomes_list = os.listdir(path)
            x3 = True
            output_file = input('Select the name of the table file with the results of the analysis: ')
            print('\n')
        else:
            print('Path to folder does not exist')
            print('Try again\n')

    for i in genomes_list:
        try:
            new_path = os.path.join(path, i)
            genome=open_file(new_path)
            error_format = False
        except:
            print('File', i, "doesn't processed! Unreadable file \n")
            error_format = True
            error_number += 1

        if error_format == False:
            forma = file_format_def(genome)
            if forma == 'error':
                print('File', i, "doesn't processed! Wrong file format")
            if forma == 'gbk':
                try:
                    fna_path = os.path.join(path, 'genomes_fna_translates')
                    if not os.path.isdir(fna_path):
                        os.mkdir(fna_path)
                    gbk_to_fna(genome, i, fna_path)
                    newfile = os.path.join(fna_path, i+'.fna')
                    print('File', newfile, 'was made!')
                    new_row = fna_in_nucleic_counter(newfile)
                    start_frame = pd.concat([start_frame, new_row], axis=1)
                    print('File', newfile, 'was processed!\n')
                    genome_number += 1
                except:
                    print('File', newfile, "doesn't processed! Inner mistake\n")
                    error_number += 1
            if forma == 'fasta':
                try:
                    new_row = fna_in_nucleic_counter(i)
                    start_frame = pd.concat([start_frame, new_row], axis=1)
                    print('File', newfile, 'was processed!\n')
                    genome_number += 1
                except:
                    print('File', newfile, "doesn't processed! Inner mistake\n")
                    error_number += 1
    if genome_number == 0:
        print('Error, no files found')
    if genome_number >= 0:
        print('\nProcessed', genome_number, 'files')
        print('\nThere were', error_number, 'errors during parsing')
        print('\nTotal processed', all_AA_count, 'sequences')
        print('\nTotal processed', all_codon_count, 'codons')
        start_frame.to_excel(output_file+'.xlsx')
        print('\nRecorded file '+output_file+'.xlsx')
        print('\n')
