import numpy as np
import math
import matplotlib.pyplot as plt
import os

def write_to_file(file_loc, data_out, separator = ' ', headers=None, append=False):
    '''inputs: file_loc-location to which to write to, data_to_write-this is a 2d array that will be written to a file
       outputs: none
       This function writes a 2d array to a file'''
       
    if not append:
        out = open(file_loc ,'w')
    else:
        out = open(file_loc ,'a+')
    
    if type(data_out) is str:
        out.write(data_out+"\n")
        
    else:
        for i  in range(len(data_out)):
            if isinstance(data_out[i], float) or type(data_out[i]) is str:
                out.write(str(data_out[i])+separator)
            else:
                for j in range(len(data_out[i])):
                    out.write(str(data_out[i][j])+separator)
                
            out.write("\n")
    out.close()
            
            
def read_file(file_loc, comments ='#', separators = [], skip_lines = [-1], header_row = None, header_marker=None):
    '''inputs: file_loc-string of file location to read in
       outputs: data-2d array of data read in and passed out
       This function reads in a dataset in a file located at file_loc and returns a 2d array of this'''
    headers = ''
    row_num = 0
    it = 0
    skipped_lines = []
    with open(file_loc) as f:
        lines = []
        for line in f:
            if row_num not in skip_lines:
                if len(separators) > 0:
                    for separator in separators:
                        line = line.replace(separator, ' ')
                if len(line.split()) > 0 and line.split()[0] != comments:
                    if header_row is not None and it == 0:
                        headers = line.rstrip()
                        
                    else:
                        lines.append(line.rstrip().split())
                    it += 1
            elif row_num in skip_lines:
                skipped_lines.append(line)
            row_num += 1
        
    if len(separators) > 0:
        for separator in separators:
            headers = headers.replace(separator, ' ')
    
    if len(headers) > 0:
        headers = headers.split()
    
    dict = {}
    if len(headers) > 0:
        for header in headers:
            dict[header] = []
    
    for i in range(len(lines)):
        for j in range(len(lines[i])):
            if len(headers) < j+1:
                if j not in dict.keys():
                    dict[j] = []
                dict[j].append(lines[i][j])
            else:
                if headers[j] not in dict.keys():
                    dict[headers[j]] == []
                dict[headers[j]].append(lines[i][j])
        
    return {"lines": lines, "dict": dict, "skipped_lines":skipped_lines}
