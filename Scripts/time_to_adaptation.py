import matplotlib as mpl
mpl.use('pdf')
import sys
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
    
plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

num = 1

x = [0.00005,0.0005,0.005,0.05,0.5]
x2 = ['0.01','0.1','1.0','10.0','100.0']
color_shape = ['go','b^','r*']

timeToAdaptation = []
numShared = []
totalNumFWAA = []
numShared_DivTotal = []
corFWAA = []
corEffect = []

for j in [1]:

	TTA = []
	num = str(j)
	recipes = ['9_'+num+'_0','9_'+num+'_1','9_'+num+'_2','9_'+num+'_3','9_'+num+'_4']
	
	for i in recipes:

		File = open("./Output1/MyRecipe"+i+"/Adaptation.txt","r")
		header = File.readline()
		FWAA_P2 = File.readline()
		FWAA_P3 = File.readline()

		TTA.append(int(File.readline().split()[0]))
		numShared_ = File.readline().split()[0]
		totalNumFWAA_ = File.readline().split()[0]

		numShared_DivTotal.append(float(numShared_ )/float(totalNumFWAA_))

		corFWAA.append(File.readline().split()[0])
		corEffect.append(File.readline().split()[0])

	timeToAdaptation.append(TTA)

timeToAdaptation = np.array(timeToAdaptation)	
means = np.mean(timeToAdaptation,axis=0)

f1, ax1 = plt.subplots()
for k in range(len(timeToAdaptation)):
	ax1.plot(x2[:-1],timeToAdaptation[k][:-1],color_shape[k])
ax1.plot(x2[:-1],means[:-1],"y--")
ax1.plot(x2[-1],10,"")

ax1.grid(axis = 'y')
ax1.set_ylabel('Time (Generations)')
ax1.set_xlabel('$M$')
ax1.set_title('Time to Adaptation')

ax1.set_yscale('log')
f1.subplots_adjust(left=.15, bottom=.16, right=.85, top=.90)
height = 6.00
width = 3.00
f1.set_size_inches(height, width)
f1.savefig("Time_Adapt.pdf")



