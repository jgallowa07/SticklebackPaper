import sys
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

num = "1"

Traces = []

recipes = ['9_'+num+'_0','9_'+num+'_1','9_'+num+'_2','9_'+num+'_3','9_'+num+'_4']

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

for i in range(5):

	c1 = []
	c2 = []
	c3 = []

	File = open("./Output1/MyRecipe"+recipes[i]+"/AveragePhenotypeThroughout.txt","r")
	first = File.readline()
	numSamplingPoints,interval,intro,end = (int(s) for s in first.split()) 

	newLakesIntroduced = ( intro / interval ) -1

	count = 0
	for line in File:
		if(count < newLakesIntroduced):
			li = line.split()
			c1.append(float(li[0]))
			c2.append(float(li[1]))
		else:
			li = line.split()
			c1.append(float(li[0]))
			c2.append(float(li[1]))
			c3.append(float(li[27]))

		count += 1;		

	col1.append(c1)
	col2.append(c2)
	col3.append(c3)

data_a = col1
data_b = col2
data_c = col3

ticks = ['0.01','0.1','1.0','10.0','100.0']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

fig1,ax1 = plt.subplots()
print(len(data_a))
print(np.array(range(len(data_a)))*3.0-0.6)

bpl = ax1.boxplot(data_a, positions=np.array(range(len(data_a)))*3.0-0.6, sym='', widths=0.4)
bpm = ax1.boxplot(data_b, positions=np.array(range(len(data_b)))*3.0, sym='', widths=0.4)
bpr = ax1.boxplot(data_c, positions=np.array(range(len(data_c)))*3.0+0.6, sym='', widths=0.4)
set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
set_box_color(bpm, '#2C7BB6')
set_box_color(bpr, '#5D7FF6')

ax1.plot([], c='#D7191C', label='Marine')
ax1.plot([], c='#2C7BB6', label='Old lakes')
ax1.plot([], c='#5D7FF6', label='New lakes')
ax1.set_ylabel('Mean trait value',fontsize=12)
ax1.set_xlabel('$M$',fontsize=12)
ax1.axhline(linestyle = '--',linewidth=2, color='black',label = 'Ancestral trait value')
ax1.axhline(linestyle = '--',y=10,linewidth=2, color='pink',label = 'Marine optimum')
ax1.axhline(linestyle = '--',y=-10,linewidth=2, color='#EE82EE',label = 'Freshwater optimum')
ax1.legend(bbox_to_anchor=(1.00, 1), loc=2, frameon=False,fontsize = 8)
plt.xticks(range(0, len(ticks) * 3, 3), ticks)
plt.yticks(np.arange(-12,12,2.0))
plt.xlim(-1.0, len(ticks)*2.7)
plt.ylim(-11, 11)
fig1.subplots_adjust(left=.15, bottom=.16, right=.75, top=.90)
width = 6.00
height = 3.00
fig1.set_size_inches(width,height)
fig1.savefig("Pheno_Dist.pdf")

