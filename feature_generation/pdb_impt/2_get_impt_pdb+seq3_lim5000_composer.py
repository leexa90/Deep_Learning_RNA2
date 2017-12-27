import numpy as np
import sys
'''
Code produces DATA dictionary, containiing
pdbID_chain
sequence
sequence resolved (X for unresolved regions)
distance matrix
'''
sys.path.append('/home/leexa/pymol/RNA/ENTIRE_RNA_ONLY/')
import alignment
pdb_seq_chain =  np.load('../fasta_chain_lim5000_composer.npy').item()
from get_impt_pdb_new_resln_module import get_clusters,get_cif,get_cif_entirePDB

RNA_puzzle = ['5T5A', '5TPY', '5K7C', '5DI4', '5DDO',
              '4XW7', '4QLM', '4LCK', '4L81', '4R4V',
              '4GXY', '4P9R', '3V7E', '3OWZ', '3P59',
              '3MEI'][0:0]

result2_resln_temp  = get_clusters()[0]
result2_resln  = []
# below used to manually remove test cases
for i in result2_resln_temp:
    in_test = False
    for j in i:
        if j[2].strip()[0:4] in RNA_puzzle:
            print j[2].strip()[0:4],i
            in_test= True
    if in_test is False: #not in test
        result2_resln += [i,]
                

def score(a,b):
    a = a.upper()
    b = b.upper()
    assert len(a) == len(b)
    for i in range(0,len(a)):
        if a[i] != 'X' and b[i] != 'X':
            if a[i] != b[i] :
                return False
    else:
        return True
        
def cheap_alignment(pdb_seq0,full_seq0):
    for i in range(0,len(full_seq0)-len(pdb_seq0)+1):
        if score(full_seq0[i:i+len(pdb_seq0)],pdb_seq0):
            #print full_seq0[i:i+len(pdb_seq0)]
            #print pdb_seq0,'\n',
            return pdb_seq0, full_seq0[i:i+len(pdb_seq0)]
def print_file(cord,chain,i):
    return ('ATOM   '+' '*(5-len(str(i)))+str(i)+"  C3'   "+str(cord[(chain,i)][0])+' A'+\
            ' '*(4-len(str(i)))+str(i)+'    '+' '*(8-len(str(cord[(chain,i)][1][0])))+\
            str(cord[(chain,i)][1][0])+' '*(8-len(str(cord[(chain,i)][1][1])))+\
            str(cord[(chain,i)][1][1])+' '*(8-len(str(cord[(chain,i)][1][2])))+\
            str(cord[(chain,i)][1][2])+'\n')
'''
1CX0_b -> 4PR6_B
2joo_w -> 4V51_AW
3I9E_B -> 4v6g_DB
3KIR_B -> 4V7J_AB
3OFR_B -> 4V7U_BB
3OFQ_B -> 4V7U_DB
'''
pdb_seq_chain['1jbr_d']='gcgcuccucaguacgagaggaaccggagcgc'
pdb_seq_chain['2dr8_b']='ggcccggggcgguucgauuccgcccugggccac'
pdb_seq_chain['3ova_c']='ggaaguagaugguucaaguccauuuacuuccacc'
composer = ['4PR6_B', '1DK1_B', '1EXD_B', '1FFY_T', '1I6U_C',
            '1JBR_D', '1KXK_A', '1LNG_B', '1MMS_C', '1MZP_B',
            '1U0B_A', '1UN6_E', '1VQO_9', '1WZ2_C', '1Z43_A',
            '1ZHO_B', '2DR8_B', '2HW8_B', '4V51_AW', '2PXB_B',
            '2PXL_B', '2QUS_A', '2QWY_A', '2R8S_R', '2V3C_M',
            '2VPL_B', '3ADB_C', '3AM1_B', '3D2V_A', '4V6G_DB',
            '3IAB_R', '3IQP_A', '3IWN_A', '4V7J_AB', '3LA5_A',
            '3NDB_M', '4V7U_DB', '4V7U_BB', '3OVA_C', '3PDR_A']
def main(ii,DATA):
    chain = ii[5:]   
    pdb_seq0, full_seq0, cord = get_cif(ii)
    
    #if cheap_alignment(pdb_seq0,full_seq0):
    if ii in ['4ADX_9','4ADX_8','4V8T_1']:
        None#no atoms
    elif ii.upper() == '3J0P_W':
        full_seq0 = pdb_seq0
        pdb_seq,full_seq = cheap_alignment(pdb_seq0,full_seq0)
##                    full_seq = full_seq0
##                elif ii.upper() == '3JAH_2':
##                    full_seq = 'GUCUCCGUAGUGUAGCUGGUAUCACGUUCGCCUAACACGCGAAAGGUCCUCGGUUCGAAACCGGGCGGAAACACCA'
##                    pdb_seq = pdb_seq0
    else:
        pdb_seq,full_seq = cheap_alignment(pdb_seq0,full_seq0)
    matrix = np.zeros((len(pdb_seq),len(pdb_seq)))/0 # get nan
    fpdb = open('../final/Z'+ii+'.pdb','w')
    fpdb2 = open('../final/fasta_'+ii+'.pdb','w')
    fpdb2.write('>'+ii+'\n')
    fpdb2.write(full_seq+'\n')
    fpdb3 = open('../final/Singlefasta_'+ii+'.pdb','w')
    fpdb3.write(full_seq+'\n')
    fpdb4 = open('../final/Doublefasta_'+ii+'.pdb','w')
    fpdb4.write(full_seq+'\n')
    fpdb4.write(full_seq+'\n')
    if len(pdb_seq) < 5000:
        for i in range(1,1+len(pdb_seq)):
            if (chain,i) in cord.keys():
                i_cord = cord[(chain,i)][1]
                fpdb.write(print_file(cord,chain,i))
                for j in range(i,1+len(pdb_seq)): 
                    if (chain,j) in cord.keys():
                        j_cord = cord[(chain,j)][1]
                        dist = sum((i_cord - j_cord)**2)**.5
                        matrix[i-1,j-1] = dist
                        matrix[j-1,i-1] = dist
        DATA[ii] = [pdb_seq,full_seq,matrix]
        completed = True
        fpdb.close()
    return DATA

have = []
DATA = {}
counter =0
import gc
for i in composer:
    ii = i.lower()
    DATA = main(ii,DATA)

for cluster in result2_resln:
    counter += 1
    completed = False
    for pdb in sorted(cluster):
        ii = pdb[2].strip()[0:4]+'_'+pdb[2].strip().split('|')[-1]
        ii = ii.lower().strip()
        if ii.upper() in composer:
            completed = True
            print completed
            print len(DATA)
    for pdb in sorted(cluster): 
        if completed is False and '+' not in pdb[2]:
            try:
                ii = pdb[2].strip()[0:4]+'_'+pdb[2].strip().split('|')[-1]
                ii = ii.lower().strip()
                DATA = main(ii,DATA)
                completed = True
    
            except TypeError:
                print ii, 'TYPE_ERROR' #cannot align
        ##        print pdb_seq0
        ##        print full_seq0
            except IOError:
                None#print ii#None#DATA[ii] = 'IO_ERROR' #pdb dont exist
            except ValueError:
                None#DATA[ii] = 'index_error' #no cordinates
                print ii
        ##        print pdb_seq0
        ##        print full_seq0
            
    if completed == False:
        print sorted(cluster),counter
        gc.collect()
np.save('../data_lim5000_nan_composer.npy',DATA)               
    
