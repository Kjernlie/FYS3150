import matplotlib.pyplot as plt
import seaborn as sns

infile = open('../Project4/probabilities.dat')

vals = []


lines = infile.readlines()
lines = lines[:-1]

for line in lines:
	#print line

	number = line.split()
	vals.append(float(number[0]))
	
#plt.hist(vals,bins=1000)
plt.plot(vals)
#sns.distplot(vals, kde=False, rug=True)
#sns.kdeplot(vals, shade=True)
#sns.distplot(vals, hist=False, rug=True)

plt.show()
