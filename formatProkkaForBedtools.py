# Author: Deepika Gunasekaran
# Title: Prokka gff file format for Bedtools
# Description: This program takes as input the prokka annotation file. This file 
# is then modified to appropriate format to serve as an input for bedtools

import sys

def fileModify(query):
    return_line = ''
    for i in query:
        line_list = i.split('\t')
        if len(line_list) == 9:
            line_list_pre = line_list[:3]
            line_list_suf = line_list[5:8]
            line_last = line_list[8]
            line = 'Synteny\t' + str(line_list[3]) + '\t' + str(line_list[4]) + '\t'
            for j in line_list_pre:
                line = line + str(j) + '\t'
            for k in line_list_suf:
                line = line + str(k) + '\t'
            line += line_last
            return_line += line
            #return_line += '\n'
    return return_line

prokka_clean = open(sys.argv[1], 'r')
output_prokka = open(sys.argv[2],'w')
prokka_gff = []
for line in prokka_clean:
    prokka_gff.append(line)
    
mod_prokka = fileModify(prokka_gff)
#print mod_prokka
output_prokka.write(mod_prokka)
prokka_clean.close()
output_prokka.close()



