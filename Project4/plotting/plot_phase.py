# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import seaborn as sns

file40 = open('../data/f/40/L40mc1e6t001.dat')
file60 = open('../data/f/60/L60mc1e6t0001.dat')
file100 = open('../data/f/100/L100mc1e6t0001.dat')
file140 = open('../data/f/140/L140mc2e6t0001.dat')

energy1 = []
energy2 = []
energy3 = []
energy4 = []
temp1 = []
temp2 = []
temp3 = []
temp4 = []
mag1 = []
mag2 = []
mag3 = []
mag4 = []
cv1 = []
cv2 = []
cv3 = []
cv4 = []
sus1 = []
sus2 = []
sus3 = []
sus4 = []


for line in file40:
	numbers1 = line.split()
	energy1.append(float(numbers1[1]))
	temp1.append(float(numbers1[0]))
	mag1.append(float(numbers1[4]))
	cv1.append(float(numbers1[2]))
	sus1.append(float(numbers1[3]))

for line in file60:
	numbers2 = line.split()
	energy2.append(float(numbers2[1]))
	temp2.append(float(numbers2[0]))
	mag2.append(float(numbers2[4]))
	cv2.append(float(numbers2[2]))
	sus2.append(float(numbers2[3]))

for line in file100:
	numbers3 = line.split()
	energy3.append(float(numbers3[1]))
	temp3.append(float(numbers3[0]))
	mag3.append(float(numbers3[4]))
	cv3.append(float(numbers3[2]))
	sus3.append(float(numbers3[3]))


for line in file140:
	numbers4 = line.split()
	energy4.append(float(numbers4[1]))
	temp4.append(float(numbers4[0]))
	mag4.append(float(numbers4[4]))
	cv4.append(float(numbers4[2]))
	sus4.append(float(numbers4[3]))

fig = plt.figure()
sub1 = fig.add_subplot(221)
sub1.plot(temp1,energy1, temp2, energy2, temp3, energy3, temp4, energy4)
sub1.set_ylabel('$\\langle E \\rangle$')
sub1.set_xlabel('$T$')
sub1.legend(['$L=40$','$L=60$','$L=100$','$L=140$'], loc='best')
sub2 = fig.add_subplot(222)
sub2.plot(temp1,cv1,temp2,cv2,temp3,cv3,temp4,cv4)
sub2.set_xlabel('$T$')
sub2.set_ylabel('$C_V$')
sub3 = fig.add_subplot(223)
sub3.plot(temp1,sus1,temp2,sus2,temp3,sus3,temp4,sus4)
sub3.set_xlabel('$T$')
sub3.set_ylabel('$ \\chi $')
sub4 = fig.add_subplot(224)
sub4.plot(temp1,mag1,temp2,mag2,temp3,mag3,temp4,mag4)
sub4.set_ylabel('$\\langle ||\\mathcal{M}| \\rangle$')
sub4.set_xlabel('$T$')
#plt.plot(temp,energy)
plt.show()



