#read in one set of runs across parameter values and plot the Fst boxplots for neutral and effect mutations


import sys
import matplotlib.pyplot as plt
import numpy as np

num = sys.argv[1]

recipes = ['8_'+num+'_0','8_'+num+'_1','8_'+num+'_2','8_'+num+'_3']
ticks = ['5e-5','5e-4','5e-3','5e-2']
titles = ['Neutral','Effect']

data_a = []
data_b = []
data_c = []

for l in range(2):

	FstMarineOrig = []
	FstMarineIntro = []
	FstOrigIntro = []

	for i in range(len(recipes)):
		
		MO1 = []
		MI1 = []
		OI1 = []

		File = open("../../Output1/MyRecipe"+recipes[i]+"/MeanFstThroughout"+titles[l]+".txt","r")
		first = File.readline()
		count = 0
		for line in File:
			#cut out the first 50 data points to cut tails from boxplot
			if(count < 50):
				count+=1
				continue
			else:
				n1,n2,n3 = (float(s) for s in line.split())
				MO1.append(n1)
				MI1.append(n2)
				OI1.append(n3)
		
		FstMarineOrig.append(MO1)
		FstMarineIntro.append(MI1)
		FstOrigIntro.append(OI1)

	data_a.append(FstMarineOrig)
	data_b.append(FstMarineIntro)
	data_c.append(FstOrigIntro)

f1, ax1 = plt.subplots(1,2,sharey=True,sharex=True,linewidth=.01)

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

for k in range(2):

	bpl = ax1[k].boxplot(data_a[k], positions=np.array(range(len(data_a[k])))*3.0-0.6, sym='', widths=0.4)
	bpm = ax1[k].boxplot(data_b[k], positions=np.array(range(len(data_b[k])))*3.0, sym='', widths=0.4)
	bpr = ax1[k].boxplot(data_c[k], positions=np.array(range(len(data_c[k])))*3.0+0.6, sym='', widths=0.4)
	set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
	set_box_color(bpm, '#2C7BB6')
	set_box_color(bpr, '#5D7FF6')

	# draw temporary red and blue lines and use them to create a legend
	ax1[k].plot([], c='#D7191C', label='Marine - Original')
	ax1[k].plot([], c='#2C7BB6', label='Marine - Introduced')
	ax1[k].plot([], c='#5D7FF6', label='Original - Introduced')
	ax1[k].legend(loc = 'upper left')

	ax1[k].set_title(titles[k] + " Mutations",fontsize='medium')
	ax1[k].set_xlabel('Migration Rate',fontsize='medium')
	ax1[k].grid(axis = 'y')

ax1[0].set_ylabel('Fst',fontsize='medium')
plt.xticks(range(0, len(ticks) * 3, 3), ticks)
plt.yticks(np.arange(0,1,0.025))
plt.xlim(-1.0, len(ticks)*2.5)
plt.ylim(0, 0.25)
plt.tight_layout()
plt.show()


