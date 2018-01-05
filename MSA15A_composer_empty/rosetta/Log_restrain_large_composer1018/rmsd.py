import os
import shutil
pdb = [x[:-4] for x in  os.listdir('.') if 'pdb' in x]
for i in pdb:
    for j in ['_0.5','_0.6','_0.75']:
        print i+j
        shutil.copy(i+'.pdb',i+j)
        shutil.copy('rmsd_XA.py',i+j)
        shutil.copy('alignment.py',i+j)
