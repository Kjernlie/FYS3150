import matplotlib.pyplot as plt
import seaborn as sns

infile1 = open('../data/b/ordered_t1/intermediate_20L_40000cycles.dat')
infile2 = open('../data/b/random_t1/intermediate_20L_40000cycles.dat')
infile3 = open('../data/b/ordered_t2.4/intermediate_20L_150000cycles.dat')
infile4 = open('../data/b/ordered_t2.4/intermediate_20L_150000cycles.dat')

mc1 = 40000*20*20
mc2 = 150000*20*20
accept1=[]
accept2=[]
accept3=[]
accept4=[]

for line in infile1:
	numbers1 = line.split()
	accept1.append(float(numbers1[5])*mc1)

for line in infile2: 
	numbers2 = line.split()
	accept2.append(float(numbers2[5])*mc1)

for line in infile3:
	numbers3 = line.split()
	accept3.append(float(numbers3[5])*mc2)

for line in infile4:
	numbers4 = line.split()
	accept4.append(float(numbers4[5])*mc2)

a = ((len(accept1)+1)/2)-1
b = ((len(accept3)+1)/2)-1

cycles1 = [100*(i+1) for i in range(len(accept1[a:]))]
cycles2 = [100*(i+1) for i in range(len(accept3[b:]))]

fig = plt.figure()
plt.suptitle('Ising model, $20\\times20$ lattice, $T=1.0$')
sub1 = fig.add_subplot(211)
sub1.plot(cycles1, accept1[a:], cycles1, accept2[a:])
sub2 = fig.add_subplot(212)
sub2.plot(cycles2, accept3[b:], cycles2, accept4[b:])

plt.show()
	
