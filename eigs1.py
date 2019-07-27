import numpy as np
from numpy import linalg as lg

import csv
import numpy as np
f=open('matrix1.csv')
csv_f=csv.reader(f)
b=[]
for row in csv_f:
	b.append(row)
a=np.zeros((len(b)-1,len(b)-1))
b.pop(0)
for i in range(len(b)):
	for j in range(len(b[0])-1):
		a[i][j]=float(b[i][j])
#print(a)
f.close()
#print(a)

w,v=lg.eig(a)
#print(w)
#print(v)

f1=open("eigenvalues1.txt","w+")
f2=open("eigenvectors1.txt","w+")

for i in range(len(w)):
	f1.write(str(w[i])+"\n")
for i in range(len(v)):
	for j in range(len(v)):
		f2.write(str(v[i][j])+"\n")
#print(w)
#print(v)
f1.close()
f2.close()
