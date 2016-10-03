import matplotlib.pylab as plt
import seaborn as sns
from math import sqrt


file = open("data_rho3_w5_int.txt")

numbers = []
rho  = []
eig1 = []
eig2 = []
eig3 = [] 

for line in file:
    numbers.append([float(i) for i in line.split()])
    
for row in numbers:
    rho.append(row[0])
    eig1.append(row[1]**2)
    eig2.append(row[2]**2)
    eig3.append(row[3]**2)


norm1 = sqrt(sum(eig1))
norm2 = sqrt(sum(eig2))
norm3 = sqrt(sum(eig3))
eig1_n = []
eig2_n = []
eig3_n = []
eig1_n = [float(i)/norm1 for i in eig1]
eig2_n = [float(i)/norm2 for i in eig2]
eig3_n = [float(i)/norm3 for i in eig3]


# -----------------------------------------------------------------------------
    
file_no = open("data_rho3_w5_noint.txt")

numbers_no = []
rho  = []
eig1_no = []
eig2_no = []
eig3_no = [] 

for line in file_no:
    numbers_no.append([float(i) for i in line.split()])
    
for row in numbers_no:
    rho.append(row[0])
    eig1_no.append(row[1]**2)
    eig2_no.append(row[2]**2)
    eig3_no.append(row[3]**2)


norm1_no = sqrt(sum(eig1_no))
norm2_no = sqrt(sum(eig2_no))
norm3_no = sqrt(sum(eig3_no))
eig1_n_no = []
eig2_n_no = []
eig3_n_no = []
eig1_n_no = [float(i)/norm1_no for i in eig1_no]
eig2_n_no = [float(i)/norm2_no for i in eig2_no]
eig3_n_no = [float(i)/norm3_no for i in eig3_no]

# --------------------------------------------------------------------------------


plt.figure()
line1, = plt.plot(rho,eig1_n,label='$\lambda_1$ w/ interaction')
line2, = plt.plot(rho,eig2_n,label='$\lambda_2$ w/ interaction')
line3, = plt.plot(rho,eig3_n,label='$\lambda_3$ w/ interaction')
line4, = plt.plot(rho,eig1_n_no,'b--',label='$\lambda_1$ w/ no interaction')
line5, = plt.plot(rho,eig2_n_no,'g--',label='$\lambda_2$ w/ no interaction')
line6, = plt.plot(rho,eig3_n_no,'r--',label='$\lambda_3$ w/ no interaction')
plt.legend(handles=[line1,line2,line3,line4,line5,line6], loc = 'best')
plt.axis([0,3,0.0, 0.030])
plt.xlabel(r'$\rho$')
plt.ylabel(r'$|u(\rho)|^2$')
plt.title(r'$\omega_r = 5.0$')
plt.show()

