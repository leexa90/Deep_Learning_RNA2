import pandas  as pd
import numpy as np
x = []
for line in open('result_comb','r'):
    x += [line.split(),]
x = pd.DataFrame(x,columns=['name','','','','','bal_acc','acc'])

x['acc'] = x['acc'].astype(np.float32)
x['bal_acc'] = x['bal_acc'].astype(np.float32)
x = x[x.acc >= 0.5][['name','bal_acc']]


def got(y,impt):
    if 'olution' in y or 'TPP' in y:
        return 0
    elif y[0:len(impt)] == impt:
        return 1
    else:
        return 0
result ={}
result_pd = pd.DataFrame(None,columns=['name','bal_acc','rank'])
for impt in [3,4,6,7,8,12,14,15,17,18,19]:
    result[impt] = x[x['name'].apply(lambda x : got(x,str(impt)))==1].copy(deep =True)
    result[impt] = result[impt].sort_values('bal_acc',ascending=False)
    result[impt]['rank'] = range(1, len(result[impt]) +1 )
    result[impt]['rank'] = result[impt]['rank'] *1.0 / len(result[impt]['rank'])
    result_pd= result_pd.append(result[impt])

def get_group(y,impt):
    if  impt.upper() in y.upper():
        return 1
    else:
        return 0
x = result_pd.copy()
result_grp  = {}
result_grp['Rnacomposer']=x[(x['name'].apply(lambda x : get_group(x,'oser'))==1)| (x['name'].apply(lambda x : get_group(x,'Adam'))==1)]
result_grp['Das'] = x[(x['name'].apply(lambda x : get_group(x,'Das'))==1) |(x['name'].apply(lambda x : get_group(x,'rhiju'))==1)]
#result_grp['major'] = x[x['name'].apply(lambda x : get_group(x,'major'))==1]
#result_grp['simrna'] = x[x['name'].apply(lambda x : get_group(x,'simrna'))==1]
#result_grp['3drna'] = x[x['name'].apply(lambda x : get_group(x,'3drna'))==1]
#result_grp['feng'] = x[x['name'].apply(lambda x : get_group(x,'feng'))==1]
#result_grp['xiao'] = x[x['name'].apply(lambda x : get_group(x,'xiao'))==1]
#result_grp['lee'] = x[x['name'].apply(lambda x : get_group(x,'lee'))==1]
result_grp['Chen'] = x[x['name'].apply(lambda x : get_group(x,'chen'))==1]
result_grp['Bujnicki'] = x[x['name'].apply(lambda x : get_group(x,'bujnicki'))==1]
result_grp['Ding'] = x[x['name'].apply(lambda x : get_group(x,'ding'))==1]
#result_grp['kevin'] = x[x['name'].apply(lambda x : get_group(x,'kevin'))==1]
result_grp['Dokholyan'] = x[x['name'].apply(lambda x : get_group(x,'_do'))==1]
#result_grp['rhiju'] = x[x['name'].apply(lambda x : get_group(x,'rhiju'))==1]
#result_grp['YagoubAli'] = x[x['name'].apply(lambda x : get_group(x,'YagoubAli'))==1]
#result_grp['mikolajczak'] = x[x['name'].apply(lambda x : get_group(x,'mikolajczak'))==1]
#result_grp['santalucia'] = x[x['name'].apply(lambda x : get_group(x,'santalucia'))==1]
result_grp['ZDeep_learning'] = x[x['name'].apply(lambda x : get_group(x,'deep'))==1]
import matplotlib.pyplot as plt
plt.hist(list(result_grp['ZDeep_learning']['rank'].values),bins=np.linspace(0,1,5),label='Deep Learning')
plt.legend()
plt.savefig('DeepLearning.png',dpi=300)
plt.close()
all_result = []
all_name= []
for i in sorted(result_grp):
    names = [x.split('_')[0] for x in result_grp[i].name]
    temp1,temp2  = [],[]
    for j in list(set(names)):
        temp = result_grp[i][result_grp[i]['name'].apply(lambda x : got(x,str(j)))==1].copy(deep =True)['rank'].values
        temp1 += [np.mean(temp),]
        temp = result_grp[i][result_grp[i]['name'].apply(lambda x : got(x,str(j)))==1].copy(deep =True)['bal_acc'].values
        temp2 += [np.mean(temp),]

        
    print i ,len(set(names)),np.mean(temp2),np.std(temp2)#,[int(xx*100) for xx in temp1]
    all_result += [temp1,]
all_name = ['Bujnicki', 'Chen', 'Das', 'Ding', 'Dokholyan', 'Rnacomposer','Deep_learning',]
plt.hist(all_result,bins=np.linspace(0,1,4),label=all_name)
plt.legend(loc=1,)
plt.ylim([0,10])
plt.ylabel('number of test cases')
plt.title('Rank distribution for each method')
plt.xticks(np.linspace(0,1,7),['','Best','','Average','','Poor',''])
plt.savefig('Compare_method.png',dpi=300)
#plt.show()
plt.close()
plt.hist(all_result,bins=np.linspace(0,1,6),label=all_name)
plt.legend(loc=1,)
plt.ylim([0,10])
plt.ylabel('number of test cases')
plt.title('Rank distribution for each method')
plt.xticks(np.linspace(0,1,11),['','Worst','','Poor','','Average','','Good','','Excellent',''][::-1])
plt.savefig('Compare_method2.png',dpi=300)
plt.close()
plt.hist(all_result,bins=np.linspace(0,1,11),label=all_name)
plt.legend(loc=1,)
plt.ylim([0,5])
plt.ylabel('number of test cases')
plt.title('Rank distribution for each method')
plt.xticks(np.linspace(0,1,21),\
           ['','Worst','','Poorer','','Poor','','Mediore','','Average','','Average','','Good','','Better',''\
            ,'Great','','Excellent',''][::-1])
plt.savefig('Compare_method3.png',dpi=300)
plt.show()
