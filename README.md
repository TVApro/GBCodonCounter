## GBCodonCounter
# Codon counter for easy work with annotated genomes

It's easy. 

You need to download one of the files and use it.

Amino_Counter is an old version for counting aminoacids in genome annotation

Codon_counter is a new version, which counting both codons and aminoacids in genome sequences and annotations 

Just install Python 3, open .py file in IDLE or another IDE,
read the instructions inside and work!


### The instruction from the Codon_Counter.py file

Analysis requires pandas, openpyxl and collections modules:

     pip3 install pandas
     pip3 install openpyxl
     pip3 install collections

The program accepts as input files:

*.gb, .gbk, .gbff and other GenBank files containing annotations

*.fna files containing nucleotide sequences of genes (be careful - the program is not yet able to adequately exclude tRNA, rRNA and ncRNA from such files)

The program generates *.fna files without RNA itself in a separate folder, they can be used later
 
### P.S.

Warning 1. The files will contain numbers with dots, which are not defined as numbers in Russian localization.

Solution 1: "Edit - Replace All" in any office program (MS Excel, OnlyOffice, LibreOffice, WPS Office etc.), replace all dots with commas

Warning 2. The error of determined percentages is set to the 6th decimal place by default.
Since the average number of codons in a bacterial genome is approximately 600,000 - 700,000, for eukaryotic comparisons
or viruses someone may need to expand or narrow the margin of error.

Solution 2: Line 279: X = format(float((X/summ)*100), '.6f')
Instead of 6, put any desired number of decimal places


Happy exploration!
