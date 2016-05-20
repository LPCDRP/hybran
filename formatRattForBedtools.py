# Author: Deepika Gunasekaran
# Title: Ratt gff file format for Bedtools
# Description: This program takes as input the ratt file which contains the 
# with no synteny with reference.This file is then modified to appropriate 
# format to serve as an input for bedtools

import sys

def fileModify(query):
    return_line = ''
    for i in query:
        line_list = i.split('\t')
        line_list_pre = line_list[:3]
        line_list_suf = line_list[5:]
        line = 'Synteny\t' + str(line_list[3]) + '\t' + str(line_list[4]) + '\t'
        for j in line_list_pre:
            line = line + str(j) + '\t'
        for k in line_list_suf:
            line = line + str(k)
        return_line += line
        #return_line += '\n'
    return return_line

ratt_clean = open(sys.argv[1], 'r')
output_ratt = open(sys.argv[2],'w')
ratt_gff = []
for line in ratt_clean:
    ratt_gff.append(line)
    
mod_ratt = fileModify(ratt_gff)
#print mod_ratt
output_ratt.write(mod_ratt)
ratt_clean.close()
output_ratt.close()