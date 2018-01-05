#import numpy as np
#data = np.load('test.npy').item()
import os
#import matplotlib.pyplot as plt
import shutil
#die #ensure that it cannot run#
val =['4v9e_aa', '5lyu_a', '4qjd_b', '4pr6_b', '5fq5_a', '4cxg_a',
                  '5m0h_a', '3amt_b', '4v8m_bd', '5x2h_b', '1e8s_c', '1c9s_w',
                  '2gtt_x', '3j0o_h', '3j45_2', '3j7r_s6', '2nz4_p', '2der_c',
                  '4cxg_2', '3p22_a', '3ivn_a', '3w3s_b', '3j0p_w', '5lzs_ii',
                  '4ug0_s6', '3d2v_a', '2csx_c', '2oiu_q', '4kzd_r', '2j28_8',
                  '5t5h_e', '1ffy_t', '5aka_7', '1pn7_c', '3j46_3', '4ue4_a',
                  '1i6u_c', '3jcs_6', '1j1u_b', '3wc1_p', '3eph_e', '2qwy_a',
                  '1un6_e', '1qzc_a', '4c4q_n', '4v6u_a1', '5xh6_b', '5mmm_z',
                  '2hw8_b', '1mj1_q', '5o60_b', '2zy6_a', '5hr6_c', '4v5z_bg',
                  '2zzm_b', '1p6v_b', '4v5z_ad', '2vpl_b', '1qzw_b', '4c7o_e',
                  '2xxa_f', '2zjr_y', '5kpy_a', '4bbl_y', '1pn8_d', '1lng_b',
                  '1m5o_b', '4kr6_d', '3nkb_b', '1gax_c', '4kr6_c', '2nue_c',
                  '4v8b_ab', '5t83_a', '3p49_a', '3izd_a', '5ktj_a', '3j9w_bb',
                  '3k0j_e', '5gap_v', '3ski_a', '2om7_g', '1ysh_b', '4v8p_b3',
                  '4aob_a', '5lzs_2', '2wwb_d', '3iab_r', '4qjh_b', '4yco_d',
                  '4tue_qv', '4kr7_x', '4adx_8', '2go5_9', '4v8m_be', '1emi_b',
                  '3jb9_c', '5e54_a', '4p5j_a', '1zc8_h', '1y26_x', '1zc8_a',
                  '1hc8_c', '3iyq_a', '5it9_i', '4wj3_q', '3suh_x', '1xjr_a',
                  '4frg_b', '1zn1_c']
val = ['5tpy_a', '5di4_a', '4p95_a', '5ddo_a', '4l81_a', '4gxy_a', '4qln_a', '3oww_a', '4xwf_a', '5ddo_b', '5t5a_a', '4r4p_a', '5k7c_a', '3v7e_c', '5dqk_a', '4lck_c', '4lck_b']
for ii in [x[:-4] for x in  os.listdir('.') if 'png' in x]:
    files = [x for x in os.listdir('.') if ii in x]
    for weights in (0.5,0.6,0.75):
        i = ii + '_' + str(weights)
        if not os.path.exists(i):
            os.makedirs(i)
        else:
            shutil.rmtree(i)
            os.makedirs(i)
        for f in files:
            if 'weights' not in f:
                if 'secstruct' in f:
                    shutil.copy(f,'%s/%s' %(i,'PA.secstruct'))
                elif 'fasta' in f:
                    shutil.copy(f,'%s/%s' %(i,'PA.fasta'))
                if 'png' in f:
                    shutil.copy(f,'%s/%s' %(i,f))
                    
            else:
                if str(weights) in f:
                    shutil.copy(f,'%s/%s' %(i,f))
                    shutil.copy(f,'%s/%s' %(i,'constraints'))
                shutil.copy('rosetta.sh','%s/%s' %(i,'rosetta.sh'))
                shutil.copy('1_automate_rna_helix.py','%s/%s' %(i,'1_automate_rna_helix.py'))
        











die
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
             
