import sys
import matplotlib.pyplot as plt
import numpy as np


num = sys.argv[1]
title = sys.argv[2]

Traces = []

recipes = ['8_'+num+'_0','8_'+num+'_1','8_'+num+'_2','8_'+num+'_3']
mig_rates = ['mig = 5e-5 (1/10 Ind/Gen)', 'mig = .0005 (1 Ind/Gen)','mig = .005 (10 Ind/Gen)','mig = .05 (100 Ind/Gen)']

xAxis = []
xAxis2 = []
col1 = []
col2 = []
col3 = []

numSamplingPoints = 0
interval = 0 
intro = 0 
end = 0
newLakesIntroduced = 0 

for i in range(4):

	c1 = []
	c2 = []
	c3 = []

	File = open("../../Output1/MyRecipe"+recipes[i]+"/AveragePhenotypeThroughout.txt","r")
	first = File.readline()
	numSamplingPoints,interval,intro,end = (int(s) for s in first.split()) 

	newLakesIntroduced = ( intro / interval ) -1

	count = 0
	for line in File:
		if(count < newLakesIntroduced):
			#n1, n2= (float(s) for s in line.split())
			li = line.split()
			c1.append(float(li[0]))
			c2.append(float(li[1]))
		else:
			li = line.split()
			c1.append(float(li[0]))
			c2.append(float(li[1]))
			c3.append(float(li[12]))

		count += 1;		

	col1.append(c1)
	col2.append(c2)
	col3.append(c3)


#xAxis = xAxis[200:1600] + xAxis[1800:3200] + xAxis[3400:4800] + xAxis[5000:]
#xAxis2 = xAxis2[100:800] + xAxis2[900:1600] + xAxis2[1700:2400] + xAxis2[2500:]
#col1 = col1[200:1600] + col1[1800:3200] + col1[3400:4800] + col1[5000:]
#col2 = col2[200:1600] + col2[1800:3200] + col2[3400:4800] + col2[5000:]
#col3 = col3[100:800] + col3[900:1600] + col3[1700:2400] + col3[2500:]

data_a = col1
data_b = col2
data_c = col3

ticks = ['5e-5','5e-4','5e-3','5e-2']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

#plt.figure()
fig1,ax1 = plt.subplots()
print(len(data_a))
print(np.array(range(len(data_a)))*3.0-0.6)

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
ax1.legend(loc = 'upper right')
ax1.set_title('Mean Phenotype Throughout Simulation')
ax1.set_ylabel('Mean Phenotype',fontsize='medium')
ax1.set_xlabel('Migration Rate',fontsize='medium')
ax1.axhline(linestyle = '--',linewidth=2, color='black',label = 'Ancestral Phenotype')
ax1.axhline(linestyle = '--',y=10,linewidth=2, color='pink',label = 'Marine Optimum')
ax1.axhline(linestyle = '--',y=-10,linewidth=2, color='#EE82EE',label = 'Freshwater Optimum')
#ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.xticks(range(0, len(ticks) * 3, 3), ticks)
plt.yticks(np.arange(-11,11,1.0))
plt.xlim(-1.0, len(ticks)*2.5)
plt.ylim(-11, 11)
plt.tight_layout()
plt.grid(axis = 'y')
plt.show()
#plt.savefig('tester.png')

