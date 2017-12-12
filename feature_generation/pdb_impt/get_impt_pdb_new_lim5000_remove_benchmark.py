result,result2,result3 =[],[],[]
RNA_puzzle = ['5T5A', '5TPY', '5K7C', '5DI4', '5DDO',
              '4XW7', '4QLM', '4LCK', '4L81', '4R4V',
              '4GXY', '4P9R', '3V7E', '3OWZ', '3P59',
              '3MEI']

for i in open('out','r'):
    #print i
    if 'ZZZ' in  i:
        if i[-1] == '\n':
            i = i[0:-1]
        result += [ int(i.split('ZZZ')[1].split('(')[0]),]
        if result[-1] > 35:
            if result[-1] < 5000:
                #print i.split('|')[0][-4:],i.split('Chain(s): ')[1].split(';')[0].strip()
                result2 += [i.split('ZZZ')[1].split(')')[1].split(','),]

temp = {}
ALL = []
for idx in range(0,len(result2)):
    i = result2[idx]
    in_test = False
    for j in i:
        if j.strip()[0:4] in RNA_puzzle:
            print j,i
            ALL += i
            in_test = True
            result2[idx] = []
            i = result2[idx]
            break
    if in_test is False:
        for j in i:
            pdb, chain = j.strip().split('|')[0].lower(),j.strip().split('|')[2]
            if j[1:5] not in temp:
                temp[j[1:5]] = [pdb+'_'+chain.lower(),]
            else:
                temp[j[1:5]] += [pdb+'_'+chain.lower(),]
print ALL      
import os
from shutil import copy2
for counter in range(0,len(result2)):
    if  sorted(result2)[counter] != []:
        i = sorted(result2)[counter]
        pdb, chain = i[0].strip().split('|')[0].lower(),i[0].strip().split('|')[2]
    ##    if '+' not in chain:
    ##        try:
    ##            copy2('./pdb/%s_%s.pdb' %(pdb,chain),'pdb_impt')
    ##        except IOError:
    ##            print pdb,pdb,pdb,pdb
    ##    else:
    ##        print i

import urllib2
def clean_string (str):
    temp = ''
    for i in str.split('|PDBID|CHAIN|SEQUENCE\n')[1]:
        if i != '\n':
            temp += i
    return temp
list_seq = [x for x in os.listdir('seq') ]
all_dict,all_dict2 = {}, {}
for pdb in sorted(temp):
    result3 += [pdb,]
    if pdb+'.fa' not in list_seq:
        print pdb
        http = 'http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=%s&compressionType=uncompressed' %pdb
        response = urllib2.urlopen(http)
        html = response.read()
        f1 = open('./seq/'+pdb+'.fa','w')
        f1.write(html)
        f1.close()
        all_dict[pdb] = html
    else:
        html = ''
        for line in open('./seq/'+pdb+'.fa','r'):
            html += line
        all_dict[pdb] = html
        got = True
        for k in temp[pdb]:
            chain = k[5:]
            for j in range(1,len(html.split('>'))):
                if pdb.upper()+':'+chain.upper() == html.split('>')[j].split('|')[0].upper():
                    all_dict2[pdb+'_'+chain] =  clean_string(html.split('>')[j])
                    got = False
        if got :
            None#print pdb,temp[pdb],html
        
import numpy as np
die
np.save('../fasta_chain_lim5000_excludeTest.npy',all_dict2)
result3 = list (set(result3))

