import matplotlib.pyplot as plt
import numpy as np
import sys

num = sys.argv[1]

Traces = []

recipes = ['8_'+num+'_0','8_'+num+'_1','8_'+num+'_2','8_'+num+'_3']
mig_rates = ['mig = 5e-5 (1/10 Ind/Gen)', 'mig = .0005 (1 Ind/Gen)','mig = .005 (10 Ind/Gen)','mig = .05 (100 Ind/Gen)']

xAxis = []
AvgFreshAllelesPerMarineInd = []
AvgFreshAllelesPerOrigInd = []
AvgFreshAllelesPerIntroInd = []	

for i in range(4):

	M = []
	O = []
	I = []

	File = open("../../Output1/MyRecipe"+recipes[i]+"/AvgFWAA.txt","r")
	first = File.readline()
	count = 0

	for line in File:
		n1,n2,n3 = (float(s) for s in line.split())
		M.append(n1)
		O.append(n2)
		I.append(n3)

	AvgFreshAllelesPerMarineInd.append(M)
	AvgFreshAllelesPerOrigInd.append(O)
	AvgFreshAllelesPerIntroInd.append(I)

data_a = AvgFreshAllelesPerMarineInd
data_b = AvgFreshAllelesPerOrigInd
data_c = AvgFreshAllelesPerIntroInd

ticks = ['5e-5','5e-4','5e-3','5e-2']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

#plt.figure()
fig1,ax1 = plt.subplots()
bpl = ax1.boxplot(data_a, positions=np.array(range(len(data_a)))*3.0-0.6, sym='', widths=0.4)
bpm = ax1.boxplot(data_b, positions=np.array(range(len(data_b)))*3.0, sym='', widths=0.4)
bpr = ax1.boxplot(data_c, positions=np.array(range(len(data_c)))*3.0+0.6, sym='', widths=0.4)
set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
set_box_color(bpm, '#2C7BB6')
set_box_color(bpr, '#5D7FF6')

# draw temporary red and blue lines and use them to create a legend
ax1.plot([], c='#D7191C', label='Marine')
ax1.plot([], c='#2C7BB6', label='Original')
ax1.plot([], c='#5D7FF6', label='Introduced')
ax1.legend(loc = 'upper left')

ax1.set_title('Mean FAA per Individual Throughout Simulation')
ax1.set_ylabel('# of FAA',fontsize='medium')
ax1.set_xlabel('Migration Rate',fontsize='medium')
plt.xticks(range(0, len(ticks) * 3, 3), ticks)
plt.yticks(np.arange(0,10,1))
plt.xlim(-1.0, len(ticks)*2.5)
plt.ylim(0, 10)
plt.tight_layout()
plt.grid(axis = 'y')
plt.show()
#plt.savefig('tester.png')
