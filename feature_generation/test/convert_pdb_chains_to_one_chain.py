f1 = open('5k7c.pdb','r')
f2 = open('5k7c_2.pdb','w')
import  numpy as np
chain='A'
counter = 1
for line in f1:
    if 'ATOM' in line and line[17:20].strip() in ['A','U','G','C'] and line[12:16].strip() == 'P'\
    and line[21].upper() == chain.upper():
        print line,
        f2.write( line[:23]+' '*(3-len(str(counter)))+str(counter)+line[26:])
        counter += 1
        if counter == 49:counter += 4
f2.close()
        
f1 = open('5di4.pdb','r')
f2 = open('5di4_2.pdb','w')
import  numpy as np
chain='A'
counter = 1
for line in f1:
    if 'ATOM' in line and line[17:20].strip() in ['A','U','G','C'] and line[12:16].strip() == 'P'\
    and line[21].upper() == chain.upper():
        print line,
        f2.write( line[:23]+' '*(3-len(str(counter)))+str(counter)+line[26:])
        counter += 1
        if counter == 52:counter += 1
        if counter == 56:counter += 1
f2.close()
