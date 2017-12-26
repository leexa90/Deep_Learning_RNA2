import numpy as np        
                
#ans  =[[2, 22], [3, 21], [4, 7], [9, 20], [10, 19], [11, 18], [14, 17]]
def recursive_solve(seq,ss,bp1, bp2,ans=[]):
    def mid(seq,ss,i): #returns areas where there is a hairpin 
        for i in range(i+1,len(ss)):
            if ss[i] == '.':
                None
            elif ss[i] == bp1:
                return False,i
            elif ss[i] == bp2:
                return True,i
    def main(seq,ss):
        if len(seq) != len(ss):
            print 'diff length'
        else:
            base_pair = []
            for i in range(0,len(seq)):
                if ss[i] == bp1:
                    if mid(seq,ss,i)[0] is True:
                        i_1 = mid(seq,ss,i)[1]
                        ss = ss[0:i]+'.'+ ss[i+1:i_1]+'.'+ss[i_1+1:len(ss)]
                        base_pair += [[i,i_1,]]
            #print ss,base_pair
            return ss,base_pair      
    if (bp1 not in ss) and (bp2 not in ss):
        #print ans
        ans = sorted(ans)
        helix = []
        j = 0
        while j in range(0,len(ans)-1):
            temp = [ans[j]]
            for i in range(j,len(ans)-1): 
                if (ans[i][0] - ans[i+1][0])**2 == 1 and (ans[i][1] - ans[i+1][1])**2 ==1 and i != len(ans) -2:
                    temp += [ans[i+1],]
                    
                elif i == len(ans) -2: # compare last pair of the list
                    #print ans[i],ans[i+1]
                    if (ans[i][0] - ans[i+1][0])**2 == 1 and (ans[i][1] - ans[i+1][1])**2 ==1 :
                        temp += [ans[i+1],]
                        helix += [temp,]
                    else: # if last element is a single basepair
                        helix += [temp,]
                        temp = [ans[i+1]]
                        helix += [temp,]
                    
                    j = i + 1
                    break                   
                else:
                    helix += [temp,]
                    j = i + 1
                    break
        if len(ans) ==1:
            helix=[ans]
        return helix
    else:
        y = main(seq,ss)
        return recursive_solve(seq,y[0],bp1, bp2,ans+y[1])
def get_dbn (seq):
    result = ''
    for i in seq:
        if i == 1:
            result += '('
        elif i == -1 :
            result += ')'
        else:
            result += '.'
    return result
def make_ss (seq):
    ans = np.zeros(shape=(len(seq),len(seq)))
    for i in recursive_solve(seq,seq,'{','}')+recursive_solve(seq,seq,'<','>')+\
        recursive_solve(seq,seq,'[',']')+recursive_solve(seq,seq,'(',')'):
        for j in i:
            ans[j[0],j[1]] = 1
            ans[j[1],j[0]] = 1
    new_seq = []
    for i in seq:
        if i in ['[','(','{','<'][1:2]:
            new_seq += [1.0,]
        elif i in [']',')','}','>'][1:2]:
            new_seq += [-1.0,]
        else : new_seq += [0.0,]
    return [ans,np.array(new_seq)]
