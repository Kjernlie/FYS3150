import matplotlib.pyplot as plt
import seaborn as sns

infile = open('../data/f/140/L140mc2e6t0001.dat')


temp=[]
energy=[]
cv=[]
sucs=[]
mag=[]

for line in infile:
	numbers = line.split()

	temp.append(float(numbers[0]))
	energy.append(float(numbers[1]))
	cv.append(float(numbers[2]))
	sucs.append(float(numbers[3]))
	mag.append(float(numbers[4]))

fig = plt.figure()
sub1 = fig.add_subplot(221)
sub1.plot(temp,energy)
sub1.set_ylabel('$\\langle E \\rangle$')
sub2 = fig.add_subplot(222)
sub2.plot(temp,cv)
sub2.set_ylabel('$\\langle C_v \\rangle$')
sub3 = fig.add_subplot(223)
sub3.plot(temp,sucs)
sub3.set_ylabel('$ \\chi $')
sub4 = fig.add_subplot(224)
sub4.plot(temp,mag)
sub4.set_ylabel('$\\langle ||\\mathcal{M}| \\rangle$')
#plt.plot(temp,energy)
plt.show()

