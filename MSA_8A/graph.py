f1= open('../MSA_8A/updates.log','r')
f2= open('../MSA_8B/updates.log','r')
f3= open('../MSA_8C/updates.log','r')
import numpy as np
x1 = [[],[]]
for line in f1:
    if '.' in line:
        temp = map(np.float, line.split())
        x1[0] += [temp[0],]
        x1[1] += [temp[1],]

x2 = [[],[]]
for line in f2:
    if '.' in line:
        temp = map(np.float, line.split())
        x2[0] += [temp[0],]
        x2[1] += [temp[1],]

x3 = [[],[]]
for line in f3:
    if '.' in line:
        temp = map(np.float, line.split())
        x3[0] += [temp[0],]
        x3[1] += [temp[1],]
import matplotlib.pyplot as plt

plt.plot(range(len(x1[0][0::2])),x1[0][0::2],'r');
plt.plot(range(len(x1[0][1::2])),x1[0][1::2],'b');

plt.plot(range(len(x2[0][0::2])),x2[0][0::2],'m');
plt.plot(range(len(x2[0][1::2])),x2[0][1::2],'c');

plt.plot(range(len(x3[0][0::2])),x3[0][0::2],'orange');
plt.plot(range(len(x3[0][1::2])),x3[0][1::2],'green');plt.show()
