#! /usr/bin/python3
# Программа для получения данных об аминокислотном составе генома или .faa файла

import re
from re import sub
import textwrap as tw
import pandas as pd
import openpyxl
import os

single_mode = False
group_mode = False
percentage = True
new_row = ()

def start():
    x1 = False
    global single_mode
    global group_mode
    global percentage
    print('\nДобро пожаловать в программу определения аминокислотного состава!')
    while x1 == False:
        q = input('Выберите режим работы: \n 1 - одиночный анализ \n 2 - групповой анализ \n 3 - выход \n\n Ответ: ')
        otvet_list = ('1', '2', '3')
        if q == '1': 
            single_mode = True
            group_mode = False
            print('Режим анализа единичного генома')
            x1 = True
        if q == '2': 
            group_mode = True
            single_mode = False
            print('Режим анализа группы геномов')
            x1 = True
        if q == '3':
            print('Выход')
            exit()
        if q not in otvet_list: 
            print('Выбрана неверная опция. Начните выбор заново\n')
    w = input('\nНеобходимо процентное отношение?\n 0 - нет\n пропуск - по умолчанию да\n\n Ответ: ')
    if w == '0':
        percentage = False
        print('Выбрано определение абсолютного количества аминокислот')
    else:
        print('Выбрано вычисление процентного отношения\n')


def open_file(seq):
    with open(seq, 'r') as file:
            text = file.read()
    print('открыт файл', file)
    return text
   
def file_format_def(text):
    forma = ''
    if text[:1] == '>':
        forma = 'fasta'
        print('Формат файла: fasta')
    if text[:5] == "LOCUS":
        forma = 'gbk'
        print('Формат файла: genbank')
    if text[:1] != '>' and text[:5] != "LOCUS":
        forma = 'error'
        print('Входной файл не соответствует формату (genbank или .faa)')
    return forma
        
def gbk_to_faa(text, file_name, path):
    # Регулярные выражения
    title_reg = re.compile(r' {5}\S+  +')
    gene_name_reg = re.compile(r'/product=\".*?\"', re.DOTALL)
    gene_id_reg = re.compile(r'/protein_id=\".*?\"', re.DOTALL)
    gene_translation_reg = re.compile(r'/translation=\"[\w\s]*\"')
    annot_split_reg = re.compile(r'/.+=\\')
    LOCUS = re.compile(r'LOCUS   .*')
    x = 0
    prot_number = 0

    def beautifull_sequence(seq):
        seq = tw.fill(seq, width=60)
        return seq

    text_list = re.split(LOCUS, text) # разбиваем геном на участки LOCUS
    for i in text_list:
        if i:
            ind_FEATURES = i.find("FEATURES ")
            ind_ORIGIN = i.find("ORIGIN")
            annotations = i[ind_FEATURES:ind_ORIGIN]
            gene_list = re.split(title_reg, annotations)
            for j in gene_list:
                    if '/product="' in j:
                        if '/translation=' in j:
                            gene_name_position = re.search(gene_name_reg, j)
                            gene_name = re.sub(r' +', ' ', (gene_name_position.group()).replace('\n', ''))
                            
                            translation_position = re.search(gene_translation_reg, j)
                            translation = re.sub(r'\n| *', '', translation_position.group())
                            prot_number += 1

                            faa_name = gene_name.replace('/product=\"', '>').replace('"', '')
                            amino_acids_draft = translation.replace('/translation=\"', '').replace('"', '')
                            amino_acids = beautifull_sequence(amino_acids_draft)
                            if gene_name == '>' or gene_name == '> ':
                                gene_name = '>Unknown_protein_' + str(lambda x: x+1)
                            
                            faafilename = file_name+'.faa'
                            new_path = os.path.join(path, faafilename)
                            if prot_number == 1:
                                with open(new_path, 'w') as faa_file:
                                    print(faa_name+'\n'+amino_acids, file=faa_file)
                            else:
                                with open(new_path, 'a') as faa_file:
                                    print(faa_name+'\n'+amino_acids, file=faa_file)
    print('Прочитано', prot_number, 'белковых последовательностей')

def aa_in_faa_counter(file):
    global percentage
    with open(file, 'r') as myfile:
        text = myfile.read()
    faa_list = {}
    ALL_AA = str()
    fasta_list = text.split(">")
    for i in fasta_list:
        if i:
            end = i.find("\n")
            fasta_name = '>'+sub('\n', '', i[:end])
            fasta_sequence = sub('\n', '', i[end:])
            faa_list[fasta_name]=fasta_sequence
            ALL_AA += fasta_sequence
    genome_name_ind = file.rfind('/')
    genome_name = file[genome_name_ind+1:]
    print('Получено', len(ALL_AA), 'аминокислотных остатков из файла', genome_name)

    Cys = ALL_AA.count('C')
    Trp = ALL_AA.count('W')
    Asp = ALL_AA.count('D')
    Phe = ALL_AA.count('F')
    Gly = ALL_AA.count('G')
    Thr = ALL_AA.count('T')
    Ser = ALL_AA.count('S')
    Met = ALL_AA.count('M')
    Ala = ALL_AA.count('A')
    Tyr = ALL_AA.count('Y')
    His = ALL_AA.count('H')
    Leu = ALL_AA.count('L')
    Glu = ALL_AA.count('E')
    Pro = ALL_AA.count('P')
    Val = ALL_AA.count('V')
    Arg = ALL_AA.count('R')
    Lys = ALL_AA.count('K')
    Asn = ALL_AA.count('N')
    Gln = ALL_AA.count('Q')
    Ile = ALL_AA.count('I')
    hidrophobic = Trp + Phe + Gly + Met + Ala + Leu + Pro + Val + Ile
    positively_charged = Lys + Arg + His
    negatively_charged = Asp + Glu
    uncharged = Gln + Asn + Tyr + Ser + Thr + Cys
    hydrophylic = Lys + Arg + His + Asp + Glu + Gln + Asn + Tyr + Ser + Thr + Cys
    summ = Cys+Trp+Asp+Phe+Gly+Thr+Ser+Met+Ala+Tyr+His+Leu+Glu+Pro+Val+Arg+Lys+Asn+Gln+Ile
    
    if percentage == False:
        new_row = pd.DataFrame({genome_name: [Cys, Trp, Asp, Phe, Gly, Thr, Ser, Met, 
                    Ala, Tyr, His, Leu, Glu, Pro, Val, Arg, Lys, Asn, Gln, 
                    Ile, '', summ, hidrophobic, positively_charged, negatively_charged, uncharged, hydrophylic]})

    if percentage == True:
        Cys_per = format(float((Cys/summ)*100), '.2f')
        Trp_per = format(float((Trp/summ)*100), '.2f')
        Asp_per = format(float((Asp/summ)*100), '.2f')
        Phe_per = format(float((Phe/summ)*100), '.2f')
        Gly_per = format(float((Gly/summ)*100), '.2f')
        Thr_per = format(float((Thr/summ)*100), '.2f')
        Ser_per = format(float((Ser/summ)*100), '.2f')
        Met_per = format(float((Met/summ)*100), '.2f')
        Ala_per = format(float((Ala/summ)*100), '.2f')
        Tyr_per = format(float((Tyr/summ)*100), '.2f')
        His_per = format(float((His/summ)*100), '.2f')
        Leu_per = format(float((Leu/summ)*100), '.2f')
        Glu_per = format(float((Glu/summ)*100), '.2f')
        Pro_per = format(float((Pro/summ)*100), '.2f')
        Val_per = format(float((Val/summ)*100), '.2f')
        Arg_per = format(float((Arg/summ)*100), '.2f')
        Lys_per = format(float((Lys/summ)*100), '.2f')
        Asn_per = format(float((Asn/summ)*100), '.2f')
        Gln_per = format(float((Gln/summ)*100), '.2f')
        Ile_per = format(float((Ile/summ)*100), '.2f')

        hidrophobic_per = format(float((hidrophobic/summ)*100), '.2f')
        positively_charged_per = format(float((positively_charged/summ)*100), '.2f')
        negatively_charged_per = format(float((negatively_charged/summ)*100), '.2f')
        uncharged_per = format(float((uncharged/summ)*100), '.2f')
        hydrophylic_per = format(float((hydrophylic/summ)*100), '.2f')

        new_row = pd.DataFrame({genome_name:  [Cys_per, Trp_per, Asp_per, Phe_per, Gly_per, Thr_per, Ser_per, Met_per, 
                    Ala_per, Tyr_per, His_per, Leu_per, Glu_per, Pro_per, Val_per, Arg_per, Lys_per, Asn_per, Gln_per, 
                    Ile_per, '', summ, hidrophobic_per, positively_charged_per, negatively_charged_per, uncharged_per, hydrophylic_per]})
    return new_row


ab = pd.DataFrame({'Имя': ['Цистеин', 'Триптофан', 'Аспартат', 'Фенилаланин', 'Глицин', 
        'Треонин', 'Серин', 'Метионин', 'Аланин', 'Тирозин', 'Гистидин', 'Лейцин',
        'Глутамат', 'Пролин', 'Валин', 'Аргинин', 'Лизин', 'Аспарагин', 'Глутамин', 'Изолейцин', '', 'Итого',
        '0', '+', '-', '+-', 'гидрофильные'], 'Сокращ.': ['Cys', 'Trp', 'Asp', 'Phe', 'Gly', 'Thr', 'Ser', 'Met',
        'Ala', 'Tyr', 'His', 'Leu', 'Glu', 'Pro', 'Val', 'Arg', 'Lys', 'Asn', 'Gln', 'Ile', '', '', '', '', '', '', ''],
        'Код': ['C', 'W', 'D', 'F', 'G', 'T', 'S', 'M', 'A', 'Y', 'H', 'L', 'E', 'P', 'V', 'R', 'K', 'N', 'Q', 'I',
        '', '', '', '', '', '', ''], 'Природа': ['+-', '0', '-', '0', '0', '+-', '+-', '0',
        '0', '+-', '+', '0', '-', '0', '0', '+', '+', '+-', '+-', '0', '', '', '', '', '', '', '']})

start()
if single_mode == True:
    genome_name_ind = file.rfind('/')
    path = file[:genome_name_ind-1]
    try:
        x2 = False
        while x2 == False:
            file = input('Введите адрес последовательности: \n')
            genome=open_file(file)
            forma = file_format_def(genome)
            if forma == 'error':
                print('\nПопробуйте ещё раз\n')
            else:
                x2 = True
        if forma == 'gbk':
            gbk_to_faa(genome, file, path)
            newfile = file+'.faa'
            new_row = aa_in_faa_counter(newfile)
        if forma == 'fasta':
            new_row = aa_in_faa_counter(file)
        df = pd.concat([ab, new_row], axis=1)
        df.to_excel('aminoacids.xlsx')
        print('\nВсё получилось, записан файл aminoacids.xlsx \n')
    except:
        print('Записать файл не удалось')


if group_mode == True:
    genomes_list = []
    genome_number = 0
    n = 0
    x3 = False
    error_format = False
    while x3 == False:
        path = input("Укажите путь к папке с геномами: ")
        if os.path.isdir(path):
            genomes_list = os.listdir(path)
            x3 = True
        else:
            print('Пути к папке не существует')
            print('Попробуйте ещё раз\n')

    for i in genomes_list:
        try:
            new_path = os.path.join(path, i)
            genome=open_file(new_path)
            
        except:
            print('Файл', i, 'не обработан!')
            print('Нечитаемый файл \n')
            error_format == True

        if error_format == False:
            forma = file_format_def(genome)
            if forma == 'error':
                print('Файл', i, 'не обработан!\n')
            if forma == 'gbk':
                try:
                    faa_path = os.path.join(path, 'genomes_faa_translates')
                    if not os.path.isdir(faa_path):
                        os.mkdir(faa_path)
                    gbk_to_faa(genome, i, faa_path)
                    newfile = os.path.join(faa_path, i+'.faa')
                    print('Файл', newfile, 'создан!')
                    new_row = aa_in_faa_counter(newfile)
                    ab = pd.concat([ab, new_row], axis=1)
                    print('Файл', newfile, 'обработан!\n')
                    n += 1
                except:
                    print('Файл', newfile, 'не обработан! Внутренняя ошибка\n')
            if forma == 'fasta':
                try:
                    new_row = aa_in_faa_counter(i)
                    ab = pd.concat([ab, new_row], axis=1)
                    print('Файл', i, 'обработан!\n')
                    n += 1
                except:
                    print('Файл', newfile, 'не обработан! Внутренняя ошибка\n')
    if n == 0:
        print('Ошибка, файлов не обнаружено')
    if n >= 0:
        print('\nОбработано', n, 'файлов')
        ab.to_excel('aminoacids.xlsx')
        print('\nВсё получилось, записан файл aminoacids.xlsx \n')
