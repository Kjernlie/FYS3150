import matplotlib.pyplot as plt
import seaborn as sns

infile = open("../Project4/intermediate_20L_250000cycles.dat")

#mc_cycles = 10000
#outs = mc_cycles/100


temp = []
energy = []
#spec_heat = []
#suscept = []
abs_mag_moment = []

for line in infile:
	numbers = line.split()	
	
	temp.append(float(numbers[0]))
	energy.append(float(numbers[1]))
#	spec_heat.append(float(numbers[2]))
#	suscept.append(float(numbers[4]))
	abs_mag_moment.append(float(numbers[4]))



a = ((len(temp)+1)/2)-1

cycles = [100*(i+1) for i in xrange(len(energy[a:]))]
fig = plt.figure()
plt.suptitle('Ising model, $20\\times20$ lattice, $T = 1.0$')
sub1 = fig.add_subplot(211)
sub1.plot(cycles, energy[a:])
#sub1.set_xlim([0,150000])
sub1.set_ylabel('$\\langle E \\rangle$')
sub2 = fig.add_subplot(212)
sub2.plot(cycles, abs_mag_moment[a:])
#sub2.set_xlim([0,150000])
sub2.set_ylabel('$\\langle |\\mathcal{M}|  \\rangle$')
sub2.set_xlabel('Number of MC cycles')
plt.show()
