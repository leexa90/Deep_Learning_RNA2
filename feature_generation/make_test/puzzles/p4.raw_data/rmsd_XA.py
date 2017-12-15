'''
The algorithm for the rotation of P unto Q starts with two sets of paired points, P and Q.
https://en.wikipedia.org/wiki/Kabsch_algorithm
'''

import sys,os
sol= [x for x in os.listdir('.') if 'solution' in x][0]
sys.path.append('../')
import alignment
import numpy as np
for num in sorted([x for x in os.listdir('.') if ('.pdb' in x and 'solution' not in x)]):
    try:
        f1= open(sol,'r')
        f2= open('%s'%num,'r')
        first_resi = True
        cord = []
        seqP = ''
        for line in f1:
            if  'ATOM' in line and line[17:20].strip() in ['A','U','G','C'] and line[12:16].strip() in ["C1*","C1'"]:
                if first_resi is True :
                    first_resi = int(line[22:26].strip()) -1
                #print line
                resNum,resiType = (int(line[22:26].strip())- first_resi,line[17:20].strip())
                X = np.array(float(line[30:38].strip()))
                Y = np.array(float(line[38:46].strip()))
                Z = np.array(float(line[46:54].strip()))
                cord += [np.array([X,Y,Z]),]
                seqP += line[17:20].strip()
        P = np.stack(cord[:])

        first_resi = True
        cord = []
        seqQ = ''
        for line in f2:
            if  'ATOM' in line and line[17:20].strip() in ['A','U','G','C'] and line[12:16].strip() in ["C1*","C1'"]:
                if first_resi is True :
                    first_resi = int(line[22:26].strip()) -1
                #print line
                resNum,resiType = (int(line[22:26].strip())- first_resi,line[17:20].strip())
                X = np.array(float(line[30:38].strip()))
                Y = np.array(float(line[38:46].strip()))
                Z = np.array(float(line[46:54].strip()))
                cord += [np.array([X,Y,Z]),]
                line2 = line
                seqQ += line[17:20].strip()
        
        Q = np.stack(cord)
        
        #remove residues not found in either
        ali = alignment.needle(seqP,seqQ)
        good_resiP,good_resiQ = [],[]
        counterP,counterQ = 0,0
        for i in range(len(ali[1])):
            if ali[0][i] == '-':
                counterQ += 1
            elif ali[1][i] == '-':
                counterP +=1
            else:
                good_resiP += [counterP,]
                counterP +=1
                good_resiQ += [counterQ,]
                counterQ += 1
        P,Q = P[good_resiP,:],Q[good_resiQ,:]
            
        
        Q = Q - np.mean(Q,0) #translating Q to P
        P = P - np.mean(P,0) #
        A = np.matmul(P.T,Q)
        V,S,W_t=np.linalg.svd(A)
        assert np.mean(A - np.matmul(np.matmul(V,np.diag(S)),W_t))  <= 0.0001
        d = np.linalg.det(np.matmul(V,W_t))
        U = np.matmul(np.matmul(V,np.diag([1,1,d])),W_t)#rotation matrix
        PU = np.matmul(P,U)
        diff = PU-Q
        rmsd = np.mean(np.sum(diff **2,1))**.5
        print num,rmsd
    except :
        print num
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
##ax.scatter3D(Q[:,0],Q[:,1],Q[:,2],c='b',linewidths=4)
##ax.scatter3D(PU[:,0],PU[:,1],PU[:,2],c='r',alpha=0.7,linewidths=4)
###ax.scatter3D(P[:,0],P[:,1],P[:,2],c='c')
##plt.show()
