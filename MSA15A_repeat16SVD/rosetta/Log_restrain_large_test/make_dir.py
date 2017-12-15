import numpy as np
data = np.load('test.npy').item()
import os
import matplotlib.pyplot as plt
for i in data:
    if i.upper() not in os.listdir('.'):
        os.makedirs(i.upper())
    data[i][0]
    f1 = open('./%s/PA.fasta' %i.upper(),'w')
    f1.write('> %s \n' %i.upper())
    f1.write(data[i][0][1].lower()+'\n')
    f1.close()
    f1 = open('./%s/PA.secstruct' %i.upper(),'w')
    f1.write(data[i][0][1].lower()+'\n')
    ss = ''
    for j in data[i][1][1]:
        if j == 0.0 :
            ss += '.'
        elif j == 1.0:
            ss += '('
        elif j == -1.0 :
            ss += ')'
    f1.write(ss +'\n')
    f1.close()
    constraints = np.argmax(data[i][-1]+np.transpose(data[i][-1],(0,2,1,3)),3)[0]
    f1 = open('./%s/constraints' %i.upper(),'w')
    f1.write('[ atompairs ]\n')
    for x in range(len(constraints)):
        for y in range(x+1,len(constraints)):
            if constraints[x,y] == 1:
                dist = '%s %s %s %s FLAT_HARMONIC %s 1 4\n' %('P',x+1,'P',y+1,12)
                f1.write( dist)
            elif constraints[x,y] == 0:
                dist = '%s %s %s %s FLAT_HARMONIC %s 1 4\n' %('P',x+1,'P',y+1,4)
                f1.write( dist)
    f, ax = plt.subplots(1,2)
    ax[0].imshow(constraints)
    a,b,c = (data[i][0][-1] < 8)*1,(data[i][0][-1] <= 15) & (data[i][0][-1] >= 8)*1,(data[i][0][-1] > 15)*1
    batch_y = np.stack((a,b,c),axis=2)
    ax[1].imshow(np.argmax(batch_y,2))
    ax[0].set_xlabel(i.upper())
    plt.savefig('./%s/%s.png' %(i.upper(),i.upper()))
    plt.clf()
    f1.close()
             
