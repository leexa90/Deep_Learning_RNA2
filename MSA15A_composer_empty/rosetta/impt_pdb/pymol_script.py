

list =['4PR6_B', '1DK1_B', '1EXD_B', '1FFY_T', '1I6U_C',
            '1JBR_D', '1KXK_A', '1LNG_B', '1MMS_C', '1MZP_B',
            '1U0B_A', '1UN6_E', '1VQO_9', '1WZ2_C', '1Z43_A',
            '1ZHO_B', '2DR8_B', '2HW8_B', '4V51_AW', '2PXB_B',
            '2PXL_B', '2QUS_A', '2QWY_A', '2R8S_R', '2V3C_M',
            '2VPL_B', '3ADB_C', '3AM1_B', '3D2V_A', '4V6G_DB',
            '3IAB_R', '3IQP_A', '3IWN_A', '4V7J_AB', '3LA5_A',
            '3NDB_M', '4V7U_DB', '4V7U_BB', '3OVA_C', '3PDR_A']


for i in list:
    name=i.split('_')[0]+i.split('_')[1]
    name2=i.split('_')[0]+'_'+i.split('_')[1]
    try:
        cmd.fetch(name)
        cmd.save(name2.lower()+'.pdb',name)
        cmd.delete(name)
        print name
    except :
        cmd.fetch(i.split('_')[0])
        cmd.select('%s and chain %s'%(i.split('_')[0],i.split('_')[1]))
        cmd.save(name2.lower()+'.pdb','sele')
        cmd.delete(i.split('_')[0])
        print i.split('_')[0]
