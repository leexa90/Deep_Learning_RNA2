import sys,os
sys.path.append('/home/leexa/pymol/RNA/ENTIRE_RNA_ONLY/')
sys.path.append('../')
import numpy as np
import alignment
import matplotlib.pyplot as plt
models = [x for x in os.listdir('.') if ('olution' not in x and 'pdb' in x)]
dictt_model = {}
for i in models[0:1]:
    dictt_model[i] = []
    f1= open(i,'r')
    seq = ''
    temp_model = {}
    for line in f1:
        if line[17:20].strip() in ['A','U','G','C']:
            if line[12:16].strip() == 'P':
                seq  += line[17:20].strip()

solution = [x for x in os.listdir('.') if ('solution' in x and 'pdb' in x)][1:]
dictt_sol = {}

def get_mat(i):
    dictt_sol[i] = []
    f1= open(i,'r')
    seq = ''
    temp_sol = {}
    counter =1
    for line in f1:
        if line[17:20].strip() in ['A','U','G','C','GUA','CYT','ADE','URI','RU3','URA','RG3','RC5']:
            if line[12:16].strip() == 'P':
                seq  += line[17:20].strip()
                cord1 = np.float(line[30:38].strip())
                cord2 = np.float(line[38:46].strip())
                cord3 = np.float(line[46:54].strip())
                temp_sol [(int(line[22:26].strip()),line[17:20].strip())]=\
                     np.array([cord1,cord2,cord3])
            if line[12:16].strip() == "C5'":
                counter += 1
    length = 126#sorted(temp_sol)[-1][0] 
    first_resi =  1
    mat = np.zeros((length,length))
    for i in range(0,len(sorted(temp_sol))):
        key_i = sorted(temp_sol)[i]
        for j in range(i+1,len(sorted(temp_sol))):
            key_j = sorted(temp_sol)[j]
            dist = (temp_sol[key_i] - temp_sol[key_j])**2
            dist = np.sum(dist)**.5
            mat[key_i[0]-1,key_j[0]-1] = dist
            mat[key_j[0]-1,key_i[0]-1] = dist
            
    data1 = mat
    a = (data1 > 16)*1
    temp_resi_map = np.stack((1-a,a),axis=2)
    #print seq
    #plt.imshow(temp_resi_map.astype(np.float32));plt.show()
    return temp_resi_map.astype(np.float32)   

result = []
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
result = []
f, ax = plt.subplots(len(models),2,figsize=(5,(len(models)*2.5)//1))
counter  = 0
for ii in sorted(models):
    answer = np.argmax(get_mat(solution[0]),2) 
    if 'major' in ii or 'das' in ii or 'chen' in ii or 'flores' in ii: # first residue does not contain P atom
        None#answer = answer[1:,1:]
    mat_model = np.argmax(get_mat(ii),2)
    score = [[],[],[]]
    for i in range(0,answer.shape[0]):
        for j in range(i+1,answer.shape[0]):
            if (np.sum(answer[i,:]) != 0 and np.sum(answer[j,:]) != 0)\
            and (np.sum(mat_model[i,:]) !=0 and np.sum(mat_model[j,:])!=0):
                if answer[i,j] == mat_model[i,j]:
                    score[answer[i,j]] += [1,]
                else:
                    score[answer[i,j]] += [0,]
    acc = np.mean([np.mean(score[0]),np.mean(score[1])])
    acc2 = np.mean(score[0]+score[1])
    print ii,answer.shape,mat_model.shape,acc,acc2
    ax[counter,0].imshow(answer)
    ax[counter,1].imshow(mat_model)
    ax[counter,1].set_xlabel(map(lambda x :str(x)[:5],[acc,acc2]))
    ax[counter,0].set_xlabel(ii)
    result += [acc,]
    counter += 1
plt.savefig('all.png');plt.close()
print np.mean(result)
print np.std(result)
result += [[0.735,.777,.831,0.870],]
plt.hist([result[:-1],result[-1]],normed=True);plt.savefig('result.png')

