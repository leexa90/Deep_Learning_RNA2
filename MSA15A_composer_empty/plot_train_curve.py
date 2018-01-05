import numpy as np
import matplotlib.pyplot as plt

f, ax = plt.subplots(1,2,figsize=(10,5));k=0

f1= open('updates.log', 'r')

temp = []
result = []
for line in f1:
    if len(line.split()) == 0:
        result += [temp,]
        temp = []
    else:
        temp += [np.float(line.split()[0]),]
ax[0].plot(range(len(result)),[x[0] for x in result],'r',label='Train')
ax[0].plot(range(len(result)),[x[1] for x in result],'b',label='Val')
ax[0].plot(range(len(result)),[x[2] for x in result],'c',label='Test')


f1= open('updates.log', 'r')

temp = []
result = []
for line in f1:
    if len(line.split()) == 0:
        result += [temp,]
        temp = []
    else:
        temp += [np.float(line.split()[1]),]
ax[1].plot(range(len(result)),[x[0] for x in result],'r',label='Train')
ax[1].plot(range(len(result)),[x[1] for x in result],'b',label='Val')
ax[1].plot(range(len(result)),[x[2] for x in result],'c',label='Test')


ax[0].legend()
ax[1].legend()
plt.savefig('train_curve.png',dpi=200)
plt.show()
