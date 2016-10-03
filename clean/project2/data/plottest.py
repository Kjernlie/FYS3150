
import matplotlib.pylab as plt
import seaborn as sns

file = open("bs.txt")

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



plt.figure()
line1, = plt.plot(rho,eig1,label='$\lambda_1$')
line2, = plt.plot(rho,eig2,label='$\lambda_2$')
line3, = plt.plot(rho,eig3,label='$\lambda_3$')
plt.legend(handles=[line1,line2,line3], loc = 'best')
#plt.axis([0,5,0.0, 0.025])
plt.xlabel(r'$\rho$')
plt.ylabel(r'Radial probability $|u(\rho)|^2$')
plt.title(r'Radial probability distributions for three lowest-lying states')
plt.show()

