import matplotlib.pyplot as plt
import seaborn as sns

infile = open("../data/b/acceptance/testing.dat")

temp = []
acceptance = []

for line in infile:
	numbers = line.split()

	temp.append(float(numbers[0]))
	acceptance.append(float(numbers[5]))


plt.plot(temp, acceptance)
plt.title('Ratio of Accepted States')
plt.xlabel('$T$')
plt.ylabel('$N_{accepted}/N_{configs}$')
plt.xlim(1.0,2.4)
plt.show()




