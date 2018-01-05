from __future__ import print_function
from __future__ import absolute_import
import warnings
import sys
sys.path.append('../../../../')
sys.path.append('../')
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import random
random.seed(1)
# something is wrong with code. Could be double refering and overwriting data1,data2etc files
data1_t = np.load('../data_lim5000_nan.npy.zip')['data_lim5000_nan.npy'].item()
data2_t = np.load('../data_lim5000_ss.npy.zip')['data_lim5000_ss'].item()
data3_t = np.load('../data_lim5000_extra.npy').item()
data4_t = np.load('../data_lim5000_MSA.npy').item()

data1 = np.load('../data_lim5000_nan_composer.npy.zip')['data_lim5000_nan_composer.npy'].item()
data2 = np.load('../data_lim5000_ss_composer.npy.zip')['data_lim5000_ss_composer'].item()
data3 = np.load('../data_lim5000_extra_composer.npy').item()
data4 = np.load('../data_lim5000_MSA_composer.npy').item()

puzzle = ['3OX0', '3OWZ', '3OXJ', '3OXE', '3OWZ', '3OWW', '3OXM', '3OWW', '3OWI', '3OXE',
          '3OXM', '3OX0', '3OXJ', '3OWI', '3OXB', '3OXD', '3OXD', '3OXB', '4LCK', '4TZV',
          '4TZZ', '4TZZ', '4TZZ', '4LCK', '4TZZ', '4TZP', '4TZW', '4TZP', '4TZZ', '4TZZ',
          '4TZV', '4TZW', '4TZZ', '4LCK', '4TZZ', '4TZP', '4LCK', '4TZP', '5EAQ', '5DQK',
          '5EAO', '5DH6', '5DI2', '5DH8', '5DH7', '5DI4', '4R4V', '5V3I', '4R4P', '3V7E',
          '3V7E', '4L81', '4OQU', '4P9R', '4P95', '4QLM', '4QLN', '4XWF', '4XW7', '4GXY',
          '5DDO', '5DDO', '5TPY','5T5A','5K7C','5DI4']

data_test = {}
data_train = {}
data_val = {}
'''
v16 - make 2 classes, <thres_distance , >=thres_distance

'''
thres_distance  = 16
def make_array(str):
    temp = [0,]*len(str)
    for i in range(len(str)):
        if str[i]=='A':
            temp[i] = [1,0,0,0,0]
        elif str[i] =='U':
            temp[i] = [0,1,0,0,0]
        elif str[i] =='G':
            temp[i] = [0,0,1,0,0]
        elif str[i] =='C':
            temp[i] = [0,0,0,1,0]
        elif str[i]=='a':
            temp[i] = [0.5,0,0,0,0]
        elif str[i] =='u':
            temp[i] = [0,0.5,0,0,0]
        elif str[i] =='g':
            temp[i] = [0,0,0.5,0,0]
        elif str[i] =='c':
            temp[i] = [0,0,0,0.5,0]
        else:
            temp[i] = [0,0,0,0,1]
    return temp
def make_array2(str):
    temp = [0,]*len(str)
    for i in range(len(str)):
        if str[i]=='*':
            temp[i] = 100
        elif str[i] == 'X':
            temp[i] = 0
        else:
            temp[i] = int(str[i])*10
    return temp
random.seed(0)
import sys
if len(sys.argv) >= 2:
    val_id = int(sys.argv[1])
else:
    val_id = 1


data1_keys_test = [ x for x in data1_t if x[0:4].upper() in puzzle]    
data1_keys = [x for x in data1 if  len(data1[x][0]) <= 500] # the short ones get val-train split
random.shuffle(data1_keys)
data1_keys_val = [x.lower() for x in ['4PR6_B', '1DK1_B', '1EXD_B', '1FFY_T', '1I6U_C',
            '1JBR_D', '1KXK_A', '1LNG_B', '1MMS_C', '1MZP_B',
            '1U0B_A', '1UN6_E', '1VQO_9', '1WZ2_C', '1Z43_A',
            '1ZHO_B', '2DR8_B', '2HW8_B', '4V51_AW', '2PXB_B',
            '2PXL_B', '2QUS_A', '2QWY_A', '2R8S_R', '2V3C_M',
            '2VPL_B', '3ADB_C', '3AM1_B', '3D2V_A', '4V6G_DB',
            '3IAB_R', '3IQP_A', '3IWN_A', '4V7J_AB', '3LA5_A',
            '3NDB_M', '4V7U_DB', '4V7U_BB', '3OVA_C', '3PDR_A']]
data1_keys_train = [x for x in data1_keys if (x not in data1_keys_val and x[0:4].upper() not in puzzle)] +\
                   [x for x in data1 if  len(data1[x][0]) > 500]
data1_keys_train = data1_keys_train[::1]
def remove_diagonals(d):
    d = d.copy()
    d[0:2,0:2] = 0
    d[-2:,-2:] = 0
    for i in range(2,len(d)-2):
        d[i-1,i] =0
        d[i,i] =0
        d[i,i-1] =0
        d[i-2,i] =0
        d[i,i] =0
        d[i,i-2] =0
        d[i+1,i] =0
        d[i,i] =0
        d[i,i+1] =0
        d[i+2,i] =0
        d[i,i] =0
        d[i,i+2] =0
    return d
    
    
print ('### TRAINSET processing ###')
for i in data1_keys_train:
    if len(data1[i][0]) >= 35:
        if i in data2.keys():
            try:
                temp1 = data1[i]
                a = (data1[i][2] > thres_distance)*1
                temp_resi_map = np.stack((a,),axis=2)
                d0= -1*(np.isnan(data1[i][2])-1) #non-nan values ==1 , nan =0
                d = remove_diagonals(d0)
                d = np.stack((d,) ,axis=2)
                d0= np.stack((d0,),axis=2)
                pair_wise_res = {('A','A') : 0, ('U','U') : 1, ('G','G') : 2, ('C','C') : 3,
                                 ('A','U') : 4, ('A','G') : 5, ('A','C') : 6,
                                 ('G','U') : 7, ('C','U') : 8,
                                 ('C','G') : 9}
                mat_pairres = np.zeros((len(temp1[1]),len(temp1[1]),10))
                for ii in range(0,len(temp1[1])):
                    for jj in range(ii,len(temp1[1])):
                        temp_paires = sorted((temp1[1][ii],temp1[1][jj]))
                        if tuple(temp_paires) in pair_wise_res:
                            index = pair_wise_res[tuple(temp_paires)]
                            mat_pairres[ii,jj,index] = 1
                            mat_pairres[jj,ii,index] = 1
                mat_pairres_con = np.zeros((len(temp1[1]),len(temp1[1]),10))
                for ii in range(0,len(temp1[1])):
                    for jj in range(ii,len(temp1[1])):
                        temp_paires = sorted((data4[i][1][ii].upper(),
                                              data4[i][1][jj].upper()))
                        if tuple(temp_paires) in pair_wise_res:
                            index = pair_wise_res[tuple(temp_paires)]
                            mat_pairres_con[ii,jj,index] = 1
                            mat_pairres_con[jj,ii,index] = 1
                data2[i][0] = np.reshape(data2[i][0],(len(data2[i][0]),len(data2[i][0]),1))
                temp2 = data2[i]
                temp2[0] = np.concatenate((data2[i][0],mat_pairres,mat_pairres_con),axis=2)
                temp3 = data3[i]
                temp4 = data4[i]
                tempF = np.concatenate((np.array(make_array(temp1[1])).T,np.array([temp2[1]]),temp3,(np.array(make_array(temp4[1])).T),np.array([make_array2(temp4[2])])))

                if len(temp2[0]) <=500:
                    data_train[i+'_ori'] = [tempF, temp1[0],temp1[1],temp_resi_map,d,temp2[1],temp2[0]]
                for window_tup in [(35,11),(50,13),(75,25),(100,33),(125,41),(150,50),(200,66),(300,100),(400,133),(500,167)]:
                    window, jump = window_tup[0], window_tup[1]
                    for repeat in range(0,len(data1[i][0]) - window+1,jump):
                        if np.mean(d0[repeat:repeat+window,repeat:repeat+window,:]) > 0.9: 
                            data_train[i+'_'+str(window)+'_'+str(repeat)] = [tempF[:,repeat:repeat+window],
                                                   temp1[0][repeat:repeat+window],
                                                   temp1[1][repeat:repeat+window],
                                                   temp_resi_map[repeat:repeat+window,repeat:repeat+window,:],
                                                   d[repeat:repeat+window,repeat:repeat+window,:],
                                                   temp2[1][repeat:repeat+window],
                                                   temp2[0][repeat:repeat+window,repeat:repeat+window]]
                if 2*len(data1[i][0])//3 <= 500:                    
                    window = 2*len(data1[i][0])//3
                    for repeat in (0,1*len(data1[i][0])//6,2*len(data1[i][0])//6):                  
                        data_train[i+'_'+str(window)+'_0.33'] = [tempF[:,repeat:repeat+window],
                                                           temp1[0][repeat:repeat+window],
                                                           temp1[1][repeat:repeat+window],
                                                           temp_resi_map[repeat:repeat+window,repeat:repeat+window,:],
                                                           d[repeat:repeat+window,repeat:repeat+window,:],
                                                           temp2[1][repeat:repeat+window],
                                                           temp2[0][repeat:repeat+window,repeat:repeat+window]]
                if 1*len(data1[i][0])//2 <= 500:
                    window = len(data1[i][0])//2
                    for repeat in (0,1*len(data1[i][0])//4,2*len(data1[i][0])//4):                  
                        data_train[i+'_'+str(window)+'_0.5'] = [tempF[:,repeat:repeat+window],
                                                           temp1[0][repeat:repeat+window],
                                                           temp1[1][repeat:repeat+window],
                                                           temp_resi_map[repeat:repeat+window,repeat:repeat+window,:],
                                                           d[repeat:repeat+window,repeat:repeat+window,:],
                                                           temp2[1][repeat:repeat+window],
                                                           temp2[0][repeat:repeat+window,repeat:repeat+window]]
                if 2*len(data1[i][0])//3 <= 500:
                    window = 3*len(data1[i][0])//4
                    for repeat in (0,1*len(data1[i][0])//8,2*len(data1[i][0])//8):                  
                        data_train[i+'_'+str(window)+'_0.75'] = [tempF[:,repeat:repeat+window],
                                                           temp1[0][repeat:repeat+window],
                                                           temp1[1][repeat:repeat+window],
                                                           temp_resi_map[repeat:repeat+window,repeat:repeat+window,:],
                                                           d[repeat:repeat+window,repeat:repeat+window,:],
                                                           temp2[1][repeat:repeat+window],
                                                           temp2[0][repeat:repeat+window,repeat:repeat+window]]
                        
            except ValueError:
                print ('%s had ValueError, probably some of inputs are of wrong dimention. Data thrown away ' %i)
                print (data2_t[i][0].shape),
print ('### DONE ###\n### VAlidation SET processing ###')
train_n = len(data_train)
for i in data1_keys_val:
    if len(data1[i][0]) >= 30:
        if i in data2.keys():
            try:
                temp1 = data1[i]
                a = (data1[i][2] > thres_distance)*1
                temp_resi_map = np.stack((a,),axis=2)
                d = -1*(np.isnan(data1[i][2])-1) #non-nan values ==1 , nan =0
                d = remove_diagonals(d) 
                d = np.stack((d,),axis=2)
                pair_wise_res = {('A','A') : 0, ('U','U') : 1, ('G','G') : 2, ('C','C') : 3,
                                 ('A','U') : 4, ('A','G') : 5, ('A','C') : 6,
                                 ('G','U') : 7, ('C','U') : 8,
                                 ('C','G') : 9}
                mat_pairres = np.zeros((len(temp1[1]),len(temp1[1]),10))
                for ii in range(0,len(temp1[1])):
                    for jj in range(ii,len(temp1[1])):
                        temp_paires = sorted((temp1[1][ii],temp1[1][jj]))
                        if tuple(temp_paires) in pair_wise_res:
                            index = pair_wise_res[tuple(temp_paires)]
                            mat_pairres[ii,jj,index] = 1
                            mat_pairres[jj,ii,index] = 1
                mat_pairres_con = np.zeros((len(temp1[1]),len(temp1[1]),10))
                for ii in range(0,len(temp1[1])):
                    for jj in range(ii,len(temp1[1])):
                        temp_paires = sorted((data4[i][1][ii].upper(),
                                              data4[i][1][jj].upper()))
                        if tuple(temp_paires) in pair_wise_res:
                            index = pair_wise_res[tuple(temp_paires)]
                            mat_pairres_con[ii,jj,index] = 1
                            mat_pairres_con[jj,ii,index] = 1
                data2[i][0] = np.reshape(data2[i][0],(len(data2[i][0]),len(data2[i][0]),1))
                temp2 = data2[i]
                temp2[0] = np.concatenate((data2[i][0],mat_pairres,mat_pairres_con),axis=2)                
                temp3 = data3[i]
                temp4 = data4[i]
                tempF = np.concatenate((np.array(make_array(temp1[1])).T,np.array([temp2[1]]),temp3,(np.array(make_array(temp4[1])).T),np.array([make_array2(temp4[2])])))
                    #         [9-features, seq, exxist_seq, cat dist_map,cat dist_map (non-zero), ss_1d, ss_2d]
                data_val[i] = [tempF, temp1[0],temp1[1],temp_resi_map,d,temp2[1],temp2[0]]
            except ValueError:
                print ('%s had ValueError, probably some of inputs are of wrong dimention. Data thrown away ' %i)
print ('### DONE ###')
print ('### BENCHMARKSET processing ###')
print ('benchmark set contains structures which were modelled in competitions, useful for benchmark') 
for i in data1_keys_test:
    if len(data1_t[i][0]) >= 30:
        if i in data2_t.keys():
            try:
                temp1 = data1_t[i]
                a = (data1_t[i][2] > thres_distance)*1
                temp_resi_map = np.stack((a,),axis=2)
                d = -1*(np.isnan(data1_t[i][2])-1) #non-nan values ==1 , nan =0
                d = remove_diagonals(d) 
                d = np.stack((d,),axis=2)
                pair_wise_res = {('A','A') : 0, ('U','U') : 1, ('G','G') : 2, ('C','C') : 3,
                                 ('A','U') : 4, ('A','G') : 5, ('A','C') : 6,
                                 ('G','U') : 7, ('C','U') : 8,
                                 ('C','G') : 9}
                mat_pairres = np.zeros((len(temp1[1]),len(temp1[1]),10))
                for ii in range(0,len(temp1[1])):
                    for jj in range(ii,len(temp1[1])):
                        temp_paires = sorted((temp1[1][ii],temp1[1][jj]))
                        if tuple(temp_paires) in pair_wise_res:
                            index = pair_wise_res[tuple(temp_paires)]
                            mat_pairres[ii,jj,index] = 1
                            mat_pairres[jj,ii,index] = 1
                mat_pairres_con = np.zeros((len(temp1[1]),len(temp1[1]),10))
                for ii in range(0,len(temp1[1])):
                    for jj in range(ii,len(temp1[1])):
                        temp_paires = sorted((data4_t[i][1][ii].upper(),
                                              data4_t[i][1][jj].upper()))
                        if tuple(temp_paires) in pair_wise_res:
                            index = pair_wise_res[tuple(temp_paires)]
                            mat_pairres_con[ii,jj,index] = 1
                            mat_pairres_con[jj,ii,index] = 1
                data2_t[i][0] = np.reshape(data2_t[i][0],(len(data2_t[i][0]),len(data2_t[i][0]),1))
                temp2 = data2_t[i]
                temp2[0] = np.concatenate((data2_t[i][0],mat_pairres,mat_pairres_con),axis=2)
                temp3 = data3_t[i]
                temp4 = data4_t[i]
                tempF = np.concatenate((np.array(make_array(temp1[1])).T,np.array([temp2[1]]),temp3,(np.array(make_array(temp4[1])).T),np.array([make_array2(temp4[2])])))
                    #         [9-features, seq, exxist_seq, cat dist_map,cat dist_map (non-zero), ss_1d, ss_2d]
                data_test[i] = [tempF, temp1[0],temp1[1],temp_resi_map,d,temp2[1],temp2[0]]
            except ValueError:
                print ('%s had ValueError, probably some of inputs are of wrong dimention. Data thrown away ' %i)

print ('### DONE ###')             
val_n = len(data_val)
test_n = len(data_test)
print ('training samples %s from %s PDBid, val samples %s from %s PDBid, test samples %s from %s PDBid,' %(train_n,len(data1_keys_train),val_n,len(data1_keys_val),test_n,len(data1_keys_test)))            
#np.save('data_all.npy',data)
dictt = {}
length = []
for i in data_train.keys():
	if tuple(i.split('_')[0:2]) not in dictt :
		dictt[tuple(i.split('_')[0:2])] = 1
	else :
		dictt[tuple(i.split('_')[0:2])] += 1
	length  += [len(data_train[i][1]),]
#plt.hist([dictt[x] for x in dictt],100);plt.show()
#plt.hist(length,100);plt.show()
# initialzie train data here
data2_x = []
data2_y = []
data2_y_nan = []
data2_y_ss = []
data2_name = []
classweight1 = 0.
classweight2 = 0.
classweight3 = 0.
for i in data_train:
        data2_x += [data_train[i][0],]
        data2_y += [data_train[i][3],]
        data2_y_nan += [data_train[i][-3],]
        data2_y_ss += [data_train[i][-1],]
        data2_name += [i,]
##        classweight1 += np.sum(data2_y[0][:,:,0] ==1)
##        classweight2 += np.sum(data2_y[0][:,:,1] ==1)
##        classweight3 += np.sum(data2_y[0][:,:,2] ==1)
##weight1 = (classweight1 / (classweight1+classweight2+classweight3))**-1
##weight2 = (classweight2 / (classweight1+classweight2+classweight3))**-1
##weight3 = (classweight3 / (classweight1+classweight2+classweight3))**-1
##weight_coeff1 = weight1 / (weight1+weight2+weight3) 
##weight_coeff2 = weight2 / (weight1+weight2+weight3) 
##weight_coeff3 = weight3 / (weight1+weight2+weight3)
##print ('Class 1 has %s percent, given %s weight' % (weight1**-1, weight_coeff1))
##print ('Class 2 has %s percent, given %s weight' % (weight2**-1, weight_coeff2))
##print ('Class 3 has %s percent, given %s weight' % (weight3**-1, weight_coeff3))
# making a batch 
data3_x = {}
data3_y = {}
data3_y_nan = {}
data3_y_ss = {}
data3_name = {}

# initialzie val data here
data2_x_val = []
data2_y_val = []
data2_y_nan_val = []
data2_y_ss_val = []
data2_name_val = []       
for i in data_val:
        data2_x_val += [data_val[i][0],]
        data2_y_val += [data_val[i][3],]
        data2_y_nan_val += [data_val[i][-3],]
        data2_y_ss_val += [data_val[i][-1],]
        data2_name_val += [i,]
print (len(data2_y),'finished intitialising total number of training/val samples')
# initialzie test data here
data2_x_test = []
data2_y_test = []
data2_y_nan_test = []
data2_y_ss_test = []
data2_name_test = []       
for i in data_test:
        data2_x_test += [data_test[i][0],]
        data2_y_test += [data_test[i][3],]
        data2_y_nan_test += [data_test[i][-3],]
        data2_y_ss_test += [data_test[i][-1],]
        data2_name_test += [i,]
    
#data2_y = np.array(data2_y)
#data2_x = np.array(data2_x)
# Create some wrappers for simplicity
print  (data2_name_val  )
epsilon = 1e-3
def batch_normalization(x,is_training=0): #normalize each channel by its statistics. 
    mean,var = tf.nn.moments(x,[0,1,2],keep_dims=False)
    scale = tf.Variable(tf.ones([x.shape[-1]]))
    beta = tf.Variable(tf.zeros([x.shape[-1]]))
    x = tf.nn.batch_normalization(x,mean,var,beta,scale,epsilon)
    return x

def conv2d(x, W, b, strides=(1,1),relu=True,padding='SAME',name=''):
    if name != '': # Conv2D wrapper, with bias and relu activation
        x = tf.nn.conv2d(x, W, strides=[1, strides[0], strides[0], 1], padding=padding , name = name)
    else: x = tf.nn.conv2d(x, W, strides=[1, strides[0], strides[0], 1], padding=padding)
    x = tf.nn.bias_add(x, b)
    x = batch_normalization(x,phase)
    #x = tf.
    if relu is True:
        return tf.nn.relu(x)
    else:
        return x
def average_pooling2d(x,window = (2,2),strides=1,padding='same'):
    x = tf.layers.average_pooling2d(x,window,strides,padding=padding)
    #x = batch_normalization (x,phase)
    return x#tf.nn.relu(x)

def max_pooling2d(x,window = (2,2),strides=1,padding='same'):
    x = tf.layers.max_pooling2d(x,window,strides,padding=padding)
    #x = batch_normalization (x,phase)
    return x#tf.nn.relu(x)
# Parameters
learning_rate = 0.0001
training_epochs = 300 
batch_size = 1
display_step = 1
n_classes =5


tf.set_random_seed(100)

# tf Graph input
x = tf.placeholder(tf.float32,[1,15,None,1],name='sequence_factors')
# pairwise distances need to define
resi_map0 = tf.placeholder(tf.float32,[1,None,None,1],name='distance_mat')
# some values are na, these will be excluded from loss
above_zero = tf.placeholder(tf.float32,[1,None,None,1],name='above_zero')
above_zero = tf.cast(above_zero,dtype=tf.float32) #TF to float
# if tf.is_nan, use convert to 0, else use original values
resi_map = tf.where(tf.is_nan(resi_map0),above_zero,resi_map0,name='sequence')
# ss_2d
ss_2d = tf.placeholder(tf.float32,[1,None,None,21],name='ss_2d')
#resi_map = tf.reshape(resi_map, shape=[1, -1, -1, 1])
#x = tf.reshape(x, shape=[1, -1, 4, 1])
rate = tf.placeholder(tf.float32,name='rate') #dropout (keep probability)
window = 5
num1 = 64/4
num2 = 96/4
num3 = 128/4
num4 = 128/4
num5 = 128/4
num6 = 128/4
num7 = 128/4
# Store layers weight & bias
weights = {
    # 1D inception layer
    '1_wc1aa': tf.Variable(tf.random_normal([15, 1, 1, num1])),
    '1_wc1ab': tf.Variable(tf.random_normal([1, window, num1, num1])),
    '1_wc1ac': tf.Variable(tf.random_normal([1, window, num1, num1])),
    '1_wc1ba': tf.Variable(tf.random_normal([15, 1, 1, num1])),
    '1_wc1bb': tf.Variable(tf.random_normal([1, window, num1, num1])),
    '1_wc1c': tf.Variable(tf.random_normal([15, 1, 1, num1])),
    '1_wc1d': tf.Variable(tf.random_normal([15, 1, 1, num1])),
    # 1D inception layer
    '2_wc1aa': tf.Variable(tf.random_normal([1, 1, num1*4, num2])),
    '2_wc1ab': tf.Variable(tf.random_normal([1, window, num2, num2])),
    '2_wc1ac': tf.Variable(tf.random_normal([1, window, num2, num2])),
    '2_wc1ba': tf.Variable(tf.random_normal([1, 1, num1*4, num2])),
    '2_wc1bb': tf.Variable(tf.random_normal([1, window, num2, num2])),
    '2_wc1c': tf.Variable(tf.random_normal([1, 1, num1*4,num2])),
    '2_wc1d': tf.Variable(tf.random_normal([1, 1, num1*4, num2])),
    '2_SS' : tf.Variable(tf.random_normal([1,1,21,43])),
    
    # 2D inception layer output 96 layer
    '3_wc1aa': tf.Variable(tf.random_normal([1, 1, num2*4+64, num3])),
    '3_wc1ab': tf.Variable(tf.random_normal([window, 1, num3, num3])),
    '3_wc1ac': tf.Variable(tf.random_normal([1, window, num3, num3])),
    '3_wc1ba': tf.Variable(tf.random_normal([1, 1, num2*4+64, num3])),
    '3_wc1bb': tf.Variable(tf.random_normal([3, 3, num3, num3])),
    '3_wc1c': tf.Variable(tf.random_normal([1, 1, num2*4+64, num3])),
    '3_wc1d': tf.Variable(tf.random_normal([1, 1, num2*4+64, num3])),


    '4_wc1aa': tf.Variable(tf.random_normal([1, 1, num3*4, num4])),
    '4_wc1ab': tf.Variable(tf.random_normal([1, window, num4, num4])),
    '4_wc1ac': tf.Variable(tf.random_normal([window, 1, num4, num4])),
    '4_wc1ba': tf.Variable(tf.random_normal([1, 1, num3*4, num4])),
    '4_wc1bb': tf.Variable(tf.random_normal([3, 3, num4, num4])),
    '4_wc1c': tf.Variable(tf.random_normal([1, 1, num3*4, num4])),
    '4_wc1d': tf.Variable(tf.random_normal([1, 1, num3*4, num4])),

    '5_wc1aa': tf.Variable(tf.random_normal([1, 1, num4*4, num5])),
    '5_wc1ab': tf.Variable(tf.random_normal([window, 1, num5, num5])),
    '5_wc1ac': tf.Variable(tf.random_normal([1, window, num5, num5])),
    '5_wc1ba': tf.Variable(tf.random_normal([1, 1, num4*4, num5])),
    '5_wc1bb': tf.Variable(tf.random_normal([3, 3, num5, num5])),
    '5_wc1c': tf.Variable(tf.random_normal([1, 1, num4*4, num5])),
    '5_wc1d': tf.Variable(tf.random_normal([1, 1, num4*4, num5])),

    '6_wc1aa': tf.Variable(tf.random_normal([1, 1, num5*4, num6])),
    '6_wc1ab': tf.Variable(tf.random_normal([1, window, num6, num6])),
    '6_wc1ac': tf.Variable(tf.random_normal([window, 1, num6, num6])),
    '6_wc1ba': tf.Variable(tf.random_normal([1, 1, num5*4, num6])),
    '6_wc1bb': tf.Variable(tf.random_normal([3, 3, num6, num6])),
    '6_wc1c': tf.Variable(tf.random_normal([1, 1, num5*4, num6])),
    '6_wc1d': tf.Variable(tf.random_normal([1, 1, num5*4, num6])),

    '7_wc1aa': tf.Variable(tf.random_normal([1, 1, num6*4, num7])),
    '7_wc1ab': tf.Variable(tf.random_normal([window, 1, num7, num7])),
    '7_wc1ac': tf.Variable(tf.random_normal([1, window, num7, num7])),
    '7_wc1ba': tf.Variable(tf.random_normal([1, 1, num6*4, num7])),
    '7_wc1bb': tf.Variable(tf.random_normal([window/2, window/2, num7, num7])),
    '7_wc1c': tf.Variable(tf.random_normal([1, 1, num6*4, num7])),
    '7_wc1d': tf.Variable(tf.random_normal([1, 1, num6*4, num7])),

    '9_out2': tf.Variable(tf.random_normal([5,5,num7*4, 1]))   
}

biases = {
    '1_bc1aa': tf.Variable(tf.random_normal([num1])),
    '1_bc1ab': tf.Variable(tf.random_normal([num1])),
    '1_bc1ac': tf.Variable(tf.random_normal([num1])),
    '1_bc1ba': tf.Variable(tf.random_normal([num1])),
    '1_bc1bb': tf.Variable(tf.random_normal([num1])),
    '1_bc1c': tf.Variable(tf.random_normal([num1])),
    '1_bc1d': tf.Variable(tf.random_normal([num1])),

    '2_bc1aa': tf.Variable(tf.random_normal([num2])),
    '2_bc1ab': tf.Variable(tf.random_normal([num2])),
    '2_bc1ac': tf.Variable(tf.random_normal([num2])),
    '2_bc1ba': tf.Variable(tf.random_normal([num2])),
    '2_bc1bb': tf.Variable(tf.random_normal([num2])),
    '2_bc1c': tf.Variable(tf.random_normal([num2])),
    '2_bc1d': tf.Variable(tf.random_normal([num2])),
    '2_SS' : tf.Variable(tf.random_normal([43])),

    '3_bc1aa': tf.Variable(tf.random_normal([num3])),
    '3_bc1ab': tf.Variable(tf.random_normal([num3])),
    '3_bc1ac': tf.Variable(tf.random_normal([num3])),
    '3_bc1ba': tf.Variable(tf.random_normal([num3])),
    '3_bc1bb': tf.Variable(tf.random_normal([num3])),
    '3_bc1c': tf.Variable(tf.random_normal([num3])),
    '3_bc1d': tf.Variable(tf.random_normal([num3])),

    '4_bc1aa': tf.Variable(tf.random_normal([num4])),
    '4_bc1ab': tf.Variable(tf.random_normal([num4])),
    '4_bc1ac': tf.Variable(tf.random_normal([num4])),
    '4_bc1ba': tf.Variable(tf.random_normal([num4])),
    '4_bc1bb': tf.Variable(tf.random_normal([num4])),
    '4_bc1c': tf.Variable(tf.random_normal([num4])),
    '4_bc1d': tf.Variable(tf.random_normal([num4])),

    '5_bc1aa': tf.Variable(tf.random_normal([num5])),
    '5_bc1ab': tf.Variable(tf.random_normal([num5])),
    '5_bc1ac': tf.Variable(tf.random_normal([num5])),
    '5_bc1ba': tf.Variable(tf.random_normal([num5])),
    '5_bc1bb': tf.Variable(tf.random_normal([num5])),
    '5_bc1c': tf.Variable(tf.random_normal([num5])),
    '5_bc1d': tf.Variable(tf.random_normal([num5])),

    '6_bc1aa': tf.Variable(tf.random_normal([num6])),
    '6_bc1ab': tf.Variable(tf.random_normal([num6])),
    '6_bc1ac': tf.Variable(tf.random_normal([num6])),
    '6_bc1ba': tf.Variable(tf.random_normal([num6])),
    '6_bc1bb': tf.Variable(tf.random_normal([num6])),
    '6_bc1c': tf.Variable(tf.random_normal([num6])),
    '6_bc1d': tf.Variable(tf.random_normal([num6])),

    '7_bc1aa': tf.Variable(tf.random_normal([num7])),
    '7_bc1ab': tf.Variable(tf.random_normal([num7])),
    '7_bc1ac': tf.Variable(tf.random_normal([num7])),
    '7_bc1ba': tf.Variable(tf.random_normal([num7])),
    '7_bc1bb': tf.Variable(tf.random_normal([num7])),
    '7_bc1c': tf.Variable(tf.random_normal([num7])),
    '7_bc1d': tf.Variable(tf.random_normal([num7])),
    
    '9_out2': tf.Variable(tf.random_normal([1]))
}

print ('number of weights')
total_parameters = 0
for variable in sorted(weights):
    # shape is an array of tf.Dimension
    shape = weights[variable].get_shape()
    variable_parameters = 1
    for dim in shape:
        variable_parameters *= dim.value
    print(variable,variable_parameters)
    total_parameters += variable_parameters
print(total_parameters, ':total parameters')

import gc
del data1,data2,data3,data4
del data1_t,data2_t,data3_t,data4_t
del data_val,data_train
gc.collect()
# construct model
# 1D first inception layer
phase = tf.placeholder(tf.bool, name='phase')
dropout = tf.placeholder(tf.float32,name='dropout')
norm_x = x#batch_normalization(x,phase)
conv1aa = conv2d(norm_x,weights['1_wc1aa'],biases['1_bc1aa'],padding='VALID')
conv1ab = conv2d(conv1aa,weights['1_wc1ab'],biases['1_bc1ab'],name='conv1ab')
conv1ac = conv2d(conv1ab,weights['1_wc1ac'],biases['1_bc1ac'],name='conv1ac')
conv1ba = conv2d(norm_x,weights['1_wc1ba'],biases['1_bc1ba'],padding='VALID')
conv1bb = conv2d(conv1ba,weights['1_wc1bb'],biases['1_bc1bb'],name='conv1bb')
conv1ca = conv2d(norm_x,weights['1_wc1c'],biases['1_bc1c'],padding='VALID')
conv1cb = average_pooling2d(conv1ca,(1,2),1,padding='same')
conv1da = conv2d(norm_x,weights['1_wc1d'],biases['1_bc1d'],padding='VALID')
conv1 = tf.concat([conv1ac,conv1bb,conv1cb,conv1da],3)
conv1 = tf.layers.dropout(conv1,rate=dropout,name='drop1',training=False)

# 1d second inception layer
conv2aa = conv2d(conv1,weights['2_wc1aa'],biases['2_bc1aa'],relu=False) #gives linear combination between channels
conv2ab = conv2d(conv2aa,weights['2_wc1ab'],biases['2_bc1ab'])
conv2ac = conv2d(conv2ab,weights['2_wc1ac'],biases['2_bc1ac'])
conv2ba = conv2d(conv1,weights['2_wc1ba'],biases['2_bc1ba'],relu=False) #gives linear combination between channels
conv2bb = conv2d(conv2ba,weights['2_wc1bb'],biases['2_bc1bb'])
conv2ca = average_pooling2d(conv1,(1,2),1,padding='same')
conv2cb = conv2d(conv1,weights['2_wc1c'],biases['2_bc1c'])
conv2da = conv2d(conv1,weights['2_wc1d'],biases['2_bc1d']) #not to over helper function with name, changed variable from conv2d to conv2da
conv2 = tf.concat([conv2ac,conv2bb,conv2cb,conv2da],3)
conv2 = tf.layers.dropout(conv2,rate=dropout,training=False)

# linear combination of 2DSS
ss_2d_2 = conv2d(ss_2d,weights['2_SS'],biases['2_SS'],padding='VALID',relu=True)
final = []
for i in range(0,0):
    mat_x = conv2[:,:,:,i]
    final += [tf.matmul(mat_x,mat_x,transpose_a=True),]
    #final += tf.reshape(final[i],(-1,100,100,1))
#final +=  [ss_2d,]
for i in range(0,num2*4): #reconstruction of image from SVD
    diag  = tf.diag(conv2[:,0,:,i])[:,:,0,:] #batch size must be 1 because of here
    ud = tf.matmul(diag[:,:,:],conv2[:,0,:,:])
    svd = tf.matmul(ud, tf.transpose(conv2[:,0,:,:],(0,2,1))) 
    final += [svd,]

y = tf.stack(final,axis=3)
y1 = tf.concat([y,ss_2d,ss_2d_2],axis = 3)
y2 = batch_normalization(y1,phase)


# 2d first inception layer
conv3aa = conv2d(y2,weights['3_wc1aa'],biases['3_bc1aa'],relu=False) #gives linear combination between channels
# 15 x 1 followed by 1 x 15 convolution == 15 x 15 convolution
conv3ab = conv2d(conv3aa,weights['3_wc1ab'],biases['3_bc1ab'],relu=False) 
conv3ac = conv2d(conv3ab,weights['3_wc1ac'],biases['3_bc1ac'],relu=False)
conv3ba = conv2d(y2,weights['3_wc1ba'],biases['3_bc1ba'],relu=False) #gives linear combination between channels
conv3bb = conv2d(conv3ba,weights['3_wc1bb'],biases['3_bc1bb'])
conv3ca = average_pooling2d(y2,(2,2),1,padding='same')
conv3cb = conv2d(y2,weights['3_wc1c'],biases['3_bc1c'])
conv3d = conv2d(y2,weights['3_wc1d'],biases['3_bc1d'])
conv3 = tf.concat([conv3ac,conv3bb,conv3cb,conv3d],3)
conv3 = tf.layers.dropout(conv3,rate=dropout,training=True)
# 2d second inception layer
conv4aa = conv2d(conv3,weights['4_wc1aa'],biases['4_bc1aa'],relu=False) #gives linear combination between channels
# 15 x 1 followed by 1 x 15 convolution == 15 x 15 convolution
conv4ab = conv2d(conv4aa,weights['4_wc1ab'],biases['4_bc1ab'],relu=False)
conv4ac = conv2d(conv4ab,weights['4_wc1ac'],biases['4_bc1ac'],relu=False)
conv4ba = conv2d(conv3,weights['4_wc1ba'],biases['4_bc1ba'],relu=False) #gives linear combination between channels
conv4bb = conv2d(conv4ba,weights['4_wc1bb'],biases['4_bc1bb'])
conv4ca = average_pooling2d(conv3,(2,2),1,padding='same')
conv4cb = conv2d(conv3,weights['4_wc1c'],biases['4_bc1c'])
conv4d = conv2d(conv3,weights['4_wc1d'],biases['4_bc1d'])
conv4 = tf.concat([conv4ac,conv4bb,conv4cb,conv4d],3)
conv4 = tf.layers.dropout(conv4,rate=dropout,training=True)
### 2d third inception layer
conv5aa = conv2d(conv4,weights['5_wc1aa'],biases['5_bc1aa'])
conv5ab = conv2d(conv5aa,weights['5_wc1ab'],biases['5_bc1ab'])
conv5ac = conv2d(conv5ab,weights['5_wc1ac'],biases['5_bc1ac'])
conv5ba = conv2d(conv4,weights['5_wc1ba'],biases['5_bc1ba'])
conv5bb = conv2d(conv5ba,weights['5_wc1bb'],biases['5_bc1bb'])
conv5ca = average_pooling2d(conv4,(2,2),1,padding='same')
conv5cb = conv2d(conv4,weights['5_wc1c'],biases['5_bc1c'])
conv5d = conv2d(conv4,weights['5_wc1d'],biases['5_bc1d'])
conv5 = tf.concat([conv5ac,conv5bb,conv5cb,conv5d],3)
conv5 = tf.layers.dropout(conv5,rate=dropout,training=True)

### 3d third inception layer
conv6aa = conv2d(conv5,weights['6_wc1aa'],biases['6_bc1aa'])
conv6ab = conv2d(conv6aa,weights['6_wc1ab'],biases['6_bc1ab'])
conv6ac = conv2d(conv6ab,weights['6_wc1ac'],biases['6_bc1ac'])
conv6ba = conv2d(conv5,weights['6_wc1ba'],biases['6_bc1ba'])
conv6bb = conv2d(conv6ba,weights['6_wc1bb'],biases['6_bc1bb'])
conv6ca = average_pooling2d(conv5,(2,2),1,padding='same')
conv6cb = conv2d(conv5,weights['6_wc1c'],biases['6_bc1c'])
conv6d = conv2d(conv5,weights['6_wc1d'],biases['6_bc1d'])
conv6 = tf.concat([conv6ac,conv6bb,conv6cb,conv6d],3)
conv6 = tf.layers.dropout(conv6,rate=dropout,training=True)

### 7d third inception layer
conv7aa = conv2d(conv6,weights['7_wc1aa'],biases['7_bc1aa'])
conv7ab = conv2d(conv7aa,weights['7_wc1ab'],biases['7_bc1ab'])
conv7ac = conv2d(conv7ab,weights['7_wc1ac'],biases['7_bc1ac'])
conv7ba = conv2d(conv6,weights['7_wc1ba'],biases['7_bc1ba'])
conv7bb = conv2d(conv7ba,weights['7_wc1bb'],biases['7_bc1bb'])
conv7ca = average_pooling2d(conv6,(2,2),1,padding='same')
conv7cb = conv2d(conv6,weights['7_wc1c'],biases['7_bc1c'])
conv7d = conv2d(conv6,weights['7_wc1d'],biases['7_bc1d'])
conv7 = tf.concat([conv7ac,conv7bb,conv7cb,conv7d],3)
conv7 = tf.layers.dropout(conv7,rate=dropout,training=True)


out = conv2d(conv7,weights['9_out2'],biases['9_out2'],relu=False)
out_softmax = tf.nn.sigmoid(out)

# Define loss and optimizer
#logit_weight = tf.constant([weight_coeff1,weight_coeff2,weight_coeff3],tf.float32)
log_loss =tf.nn.sigmoid_cross_entropy_with_logits(logits = out, labels = resi_map )
cost =  tf.reduce_mean(tf.multiply(log_loss,above_zero)) # this is masking the nan and diagonals in the loss
learning_rate = tf.Variable(0,dtype= np.float32)
update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
with tf.control_dependencies(update_ops):
    extra_optimizer = tf.train.MomentumOptimizer(learning_rate=learning_rate,momentum=0.5).minimize(cost)
def accuracy(mat_model,answer):
    mat_model = np.reshape(mat_model,(mat_model.shape[1],mat_model.shape[1]))
    answer = np.reshape(answer,(answer.shape[1],answer.shape[1]))
    score = [[0],[0],[0]]
    for i in range(0,answer.shape[1]):
        if np.sum(answer[i,:]) != 0:
            for j in range(i+1,answer.shape[1]):
                if np.sum(answer[j,:]) != 0:
                    if answer[i,j] == mat_model[i,j]:
                        score[answer[i,j]] += [1,]
                    else:
                        score[answer[i,j]] += [0,]
    return np.mean([np.mean(score[0]),np.mean(score[1])])

import os
saver = tf.train.Saver()
# Initializing the variables
init = tf.global_variables_initializer();sess = tf.Session();sess.run(init)
saved_files = [xxx[:-5] for xxx in os.listdir('.') if (xxx.startswith('model') and xxx.endswith('.ckpt.meta'))]
next_epoch = 0
if len(saved_files) >= 1:
    last_file = sorted(saved_files,key=lambda x : int(x.split('_')[-2]))[-1]
    next_epoch = int(last_file.split('_')[-2])
    print ('starting from :' ,last_file)
    saver.restore(sess,'./'+last_file)# Training cycle
    
result = {}
random.seed(0)
shuffle = range(len(data2_x))
random.shuffle(shuffle)
text = ''
tf.summary.FileWriter('logs', graph=tf.get_default_graph())
for epoch in range(next_epoch,training_epochs):
    counter = 0
    avg_cost = []
    val_cost = []
    test_cost = []
    train_acc = []
    val_acc = []
    test_acc = []
    total_batch = train_n#int(mnist.train.num_examples/batch_size)
    # Loop over all batches
    shuffle = range(len(data2_x))
    random.shuffle(shuffle)
    num = 1
    if True:
        for batch in range(0,len(shuffle),num):
            batch_list = shuffle[batch:batch+num] 
            counter += 1
            if epoch %2 == 0:
                lr = 1+np.cos(1.0*batch*3.142/len(shuffle))
            elif epoch %2 == 1:
                lr = 1+np.cos(1.0*batch*3.142/len(shuffle))
            if epoch < training_epochs//2:
                lr = lr/10
            elif epoch < 3*training_epochs//4:
                lr = lr/100
            else:
                lr = lr/1000          
            batch_x = np.array([[data2_x[i]] for i in batch_list])
            batch_y = np.array([data2_y[i]for i in batch_list])
            batch_y_nan = np.array([data2_y_nan[i]  for i in batch_list])
            batch_y_ss = np.array([data2_y_ss[i]  for i in batch_list ])
            batch_x = np.swapaxes(np.swapaxes(batch_x,1,3),1,2)
            # Run optimization op (backprop) and cost op (to get loss value)
            _, c = sess.run([extra_optimizer, cost], feed_dict={x: batch_x,
                                                          resi_map0: batch_y,
                                                          above_zero : batch_y_nan,
                                                          ss_2d : batch_y_ss,
                                                          phase : True, learning_rate : lr, dropout : 0.11111111111})

    if True:
        k,lr = 0,0
        val_acc = []
        train_acc = []
        avg_cost  = []
        for i in range(0,len(data2_x)):
            batch_x, batch_y = np.array([[data2_x[i],],]),np.array([data2_y[i],])
            batch_y_nan,batch_y_ss = np.array([data2_y_nan[i]]),np.array([data2_y_ss[i]])
            batch_x = np.swapaxes(np.swapaxes(batch_x,1,3),1,2)
            pred = sess.run( out_softmax, feed_dict={x: batch_x,resi_map0: batch_y,
                                                         above_zero : batch_y_nan, ss_2d : batch_y_ss,
                                                            phase : False, learning_rate : lr, dropout : 0})
            c  = sess.run( cost, feed_dict={x: batch_x,resi_map0: batch_y,
                                                         above_zero : batch_y_nan, ss_2d : batch_y_ss,
                                                            phase : False, learning_rate : lr, dropout : 0})
            for k in range(len(batch_y)):
                train_acc += [accuracy((pred[k]+np.transpose(pred[k],(1,0,2)))//1,batch_y[k]),]
            # Compute average loss
            avg_cost += [c,]
        for i in range(len(data2_x_val)):
                batch_x, batch_y = np.array([[data2_x_val[i],],]),np.array([data2_y_val[i],])
                batch_y_nan,batch_y_ss = np.array([data2_y_nan_val[i]]),np.array([data2_y_ss_val[i]])
                batch_x = np.swapaxes(np.swapaxes(batch_x,1,3),1,2)
                cost_i  = sess.run( cost, feed_dict={x: batch_x,resi_map0: batch_y,
                                                     above_zero : batch_y_nan, ss_2d : batch_y_ss,
                                                     phase : False, learning_rate : lr, dropout : 0})
                val_cost += [cost_i,]
                pred =sess.run( out_softmax, feed_dict={x: batch_x,resi_map0: batch_y,
                                                     above_zero : batch_y_nan, ss_2d : batch_y_ss,
                                                        phase : False, learning_rate : lr, dropout : 0})
                val_acc += [accuracy((pred[k]+np.transpose(pred[k],(1,0,2)))//1,batch_y[k]),]
                if False:
                    temp_pred = pred[k]+np.transpose(pred[k],(1,0,2))
                    f, ax = plt.subplots(1,5,figsize=(19,5));k=0
                    ax[0].imshow(temp_pred[:,:,0]>=1)
                    ax[1].imshow(temp_pred[:,:,0]>=1.6)
                    ax[2].imshow(temp_pred[:,:,0]>=1.5)
                    ax[-2].imshow(temp_pred[:,:,0] *200//20)
                    ax[-1].imshow(batch_y[k,:,:,0]>=1)
                    ax[0].set_xlabel('pred bal_acc=%s (thres-50)'%np.round(accuracy(temp_pred[:,:,0]>=1,batch_y[k,:,:,0]>=1),2))
                    ax[1].set_xlabel('pred bal_acc=%s (thres-20)'%np.round(accuracy(temp_pred[:,:,0]>=1.6,batch_y[k,:,:,0]>=1),2))
                    ax[2].set_xlabel('pred bal_acc=%s (thres-25)'%np.round(accuracy(temp_pred[:,:,0]>=1.5,batch_y[k,:,:,0]>=1),2))
                    ax[-2].set_xlabel('probabilities logloss=%s' %cost_i)
                    ax[-1].set_xlabel('actual')
                    plt.savefig(   'VAL/'+ data2_name_val[i]+'.png');plt.close()
        test_acc = []
        for i in range(len(data2_x_test)):
                batch_x, batch_y = np.array([[data2_x_test[i],],]),np.array([data2_y_test[i],])
                batch_y_nan,batch_y_ss = np.array([data2_y_nan_test[i]]),np.array([data2_y_ss_test[i]])
                batch_x = np.swapaxes(np.swapaxes(batch_x,1,3),1,2)
                cost_i  = sess.run( cost, feed_dict={x: batch_x,resi_map0: batch_y,
                                                     above_zero : batch_y_nan, ss_2d : batch_y_ss,
                                                     phase : False, learning_rate : lr, dropout : 0})
                test_cost += [cost_i,]
                pred =sess.run( out_softmax, feed_dict={x: batch_x,resi_map0: batch_y,
                                                     above_zero : batch_y_nan, ss_2d : batch_y_ss,
                                                        phase : False, learning_rate : lr, dropout : 0})
                test_acc += [accuracy((pred[k]+np.transpose(pred[k],(1,0,2)))//1,batch_y[k]),]
                if False:
                    temp_pred = pred[k]+np.transpose(pred[k],(1,0,2))
                    f, ax = plt.subplots(1,5,figsize=(19,5));k=0
                    ax[0].imshow(temp_pred[:,:,0]>=1)
                    ax[1].imshow(temp_pred[:,:,0]>=1.6)
                    ax[2].imshow(temp_pred[:,:,0]>=1.5)
                    ax[-2].imshow(temp_pred[:,:,0] *200//20)
                    ax[-1].imshow(batch_y[k,:,:,0]>=1)
                    ax[0].set_xlabel('pred bal_acc=%s (thres-50)'%np.round(accuracy(temp_pred[:,:,0]>=1,batch_y[k,:,:,0]>=1),2))
                    ax[1].set_xlabel('pred bal_acc=%s (thres-40)'%np.round(accuracy(temp_pred[:,:,0]>=1.2,batch_y[k,:,:,0]>=1),2))
                    ax[2].set_xlabel('pred bal_acc=%s (thres-25)'%np.round(accuracy(temp_pred[:,:,0]>=1.5,batch_y[k,:,:,0]>=1),2))
                    ax[-2].set_xlabel('probabilities logloss=%s' %cost_i)
                    ax[-1].set_xlabel('actual')
                    plt.savefig( data2_name_test[i]+'_%s.png'%thres_distance);plt.close()    # Display logs per epoch step
    f1 = open('updates.log','w')
    text += str(np.mean(avg_cost))+'  '+str(np.mean(train_acc))+'\n'
    text += str(np.mean(val_cost))+'  '+str(np.mean(val_acc))+'\n'
    text += str(np.mean(test_cost))+'  '+str(np.mean(test_acc))+'\n\n'
    print ("Epoch:", '%04d' % (epoch+1), "cost=", 
            "{:.9f}".format(np.mean(avg_cost)),np.mean(train_acc))
    print ("Epoch:", '%04d' % (epoch+1), "cost=", 
            "{:.9f}".format(np.mean(val_cost)),np.mean(val_acc))
    print ("Epoch:", '%04d' % (epoch+1), "cost=", 
            "{:.9f}".format(np.mean(val_cost)),np.mean(test_acc))
    f1.write(text)
    f1.close()
    save_path = saver.save(sess,'./model300_reweigh_loss_%s_%s.ckpt' %(epoch,int(100*np.mean(val_acc))))
    result[epoch] = [avg_cost,val_cost]

print ("Optimization Finished!")

plt.plot(range(0,training_epochs),[result[i][0] for i in result],label='Train')
plt.plot(range(0,training_epochs),[result[i][1] for i in result],label='Val')
plt.legend();plt.ylabel('Logloss cost') ; plt.xlabel('epoch')
plt.savefig('Train_curveLg.png')

