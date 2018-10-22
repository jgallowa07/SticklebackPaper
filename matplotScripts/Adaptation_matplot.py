

#read in the arguments from each run across paramters values and plot all the time to adaptations.

import sys
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

num = sys.argv[1]

x = [0.00005,0.0005,0.005,0.05]
x2 = ['5e-5','5e-4','5e-3','5e-2']
color_shape = ['go','b^','rO']

timeToAdaptation = []
numShared = []
totalNumFWAA = []
numShared_DivTotal = []
corFWAA = []
corEffect = []

for j in sys.argv[1:]:

	TTA = []
	num = j
	recipes = ['8_'+num+'_0','8_'+num+'_1','8_'+num+'_2','8_'+num+'_3']
	
	for i in recipes:

		File = open("../../Output1/MyRecipe"+i+"/Adaptation.txt","r")
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
	ax1.plot(x2,timeToAdaptation[k],color_shape[k])
ax1.plot(x2,means,"y--")

ax1.grid(axis = 'y')
ax1.set_ylabel('Time (Generation)')
ax1.set_xlabel('Migration Rate')
ax1.set_title('Time to Adaptation')

ax1.set_yscale('log')
plt.show()



