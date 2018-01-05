import os
import numpy as np
pdb = [x[:-4] for x in  os.listdir('.') if ('.pdb' in x[-4:] and '_' in x)]
dbn = [x[:-17] for x in  os.listdir('.') if '.dbn' in x[-4:]]
data = np.load('../../../data_lim5000_nan_composer.npy.zip')['data_lim5000_nan_composer.npy'].item()
data1 = np.load('../../../data_lim5000_ss_composer.npy.zip')['data_lim5000_ss_composer.npy'].item()
dictt = {}
for i in sorted(pdb):
    if i not in dbn:
        print i
        dictt[i] =  data1[i]
#        cmd.load(i+'.pdb')

import matplotlib.pyplot as plt

import  numpy as  np

def get_bp_dbn(x,result=[],debug=False):
    prev = 0
    got_bracket = False
    for i in range(len(x)):
        if x[i] == '(' :
            got_bracket = True
            prev = i
        elif x[prev] == '(' and x[i] == ')':
            got_bracket = True
            result += [[prev,i],]
            x = x[0:prev]+':'+x[prev+1:i]+':'+x[i+1:]
            if debug is True:
                print x
            return get_bp_dbn(x,result)
    if got_bracket == False:
        return result
                  
import os,sys

# make fasta file
f_fasta = open('../RNA_lim5000_composer_newSS.fa','w')
all_fasta = [x for x in  os.listdir('.') if ('.dbn' in x[-4:] and '_' in x)]
for ii in sorted(all_fasta):
    unimportant = ''
    f1 = open(ii,'r')
    counter = 0
    for line in f1:
        if counter <= 1:
            f_fasta.write(line)
        else:
            break
        counter += 1
f_fasta.close()


names = {}

def array_to_ss(x):
    result = ''
    for i in x:
        if i == 1:
            result += '('
        elif i == -1:
            result += ')'
        else:
            result += '.'
    return result
import alignment
def get_new_SS(oriSeq,oriSS,dssrSeq,dssrSS):
    
    
for ii in sorted(all_fasta):
    unimportant = ''
    f1 = open(ii,'r')
    for line in f1:
        line = line.upper()
        if '>' in line or 'A' in line or 'U' in line or 'G' in line or 'C' in line:
            unimportant += line
        else:
            if line[-1] == '\n':
                line= line.split()[0]
                x = line
                
            if len(x) < 5000:
                try:
                    y = np.zeros((len(x),len(x)))
                    name_y = np.zeros(len(x))
                    if len(get_bp_dbn(x,[],False)) == 0:
                        print ii,x#,unimportant
                    else:
                        print data[ii[:-17]][1]
                        print array_to_ss(data1[ii[:-17]][1])
                        print unimportant.split('\n')[1]
                        print x,'\n'
                    for i in get_bp_dbn(x,[],False):
                        y[i[0],i[1]] = 1
                        y[i[1],i[0]] = 1
                        name_y[i[0]],name_y[i[1]] = 1,-1
                    dictt[ii[:-17]] = [y,name_y]
                    if len(x) > 399:
                        None#sprint x
                        #plt.imshow(y);plt.show()
                    break #only want the first SS, no need the mean or average bla bla
                except RuntimeError:
                    print ii, 'runtime error'

np.save('../data_lim5000_ss_composer2.npy',dictt)

