import matplotlib.pyplot as plt
import seaborn as sns

infile = open("../data/b/prob/t1/system/intermediate_20L_250000cycles.dat")

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
        energy.append(float(numbers[7]))



a = ((len(temp)+1)/2)-1


sns.distplot(energy[a:],bins=200)
plt.title('Probability distribution of $\\langle E \\rangle$, $T=1.0$')
plt.ylabel('Nr. of states')
plt.xlabel('$\\langle E \\rangle$')

#print t[-1]

plt.show()
