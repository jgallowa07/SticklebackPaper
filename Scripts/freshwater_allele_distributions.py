import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10)

num = '1'

Traces = []

recipes = ['9_'+num+'_0','9_'+num+'_1','9_'+num+'_2','9_'+num+'_3','9_'+num+'_4']
mig_rates = ["$0.01$","$0.1$","$1.0$","$10.0$","$100.0$"]

def getData(percentage=False):

    filename = "/AvgFWAA"
    if(percentage):
        filename += "_divTotal"
    filename += ".txt"

    AvgFreshAllelesPerMarineInd = []
    AvgFreshAllelesPerOrigInd = []
    AvgFreshAllelesPerIntroInd = []	

    for i in range(5):

        M = []
        O = []
        I = []

        File = open("./Output1/MyRecipe"+recipes[i]+filename,"r")
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

    return AvgFreshAllelesPerMarineInd,AvgFreshAllelesPerOrigInd,AvgFreshAllelesPerIntroInd

def getDataTotal():

    totalNumFreshwaterAdaptedAlleles = []

    for i in range(5):
        T = []

        total_File = open("./Output1/MyRecipe"+recipes[i]+"/numfreshAlleles.txt","r")
        first = total_File.readline()
        header = first.split()
        numSamples = int(header[0])
        half = numSamples / 2
        count = 0
        for line in total_File:
            if count > half:
                T.append(int(line.split()[0]))
            
            count += 1	

        totalNumFreshwaterAdaptedAlleles.append(T)

    return totalNumFreshwaterAdaptedAlleles
    
    

Data1,Data2,Data3 = getData()
DataP1,DataP2,DataP3 = getData(percentage=True)
DataT = getDataTotal()

ticks = ['5e-5','5e-4','5e-3','5e-2','5e-1']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

fig,axes = plt.subplots(nrows = 1,ncols = 3,sharex=True)
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
bpl = ax1.boxplot(Data1, positions=np.array(range(len(Data1)))*3.0-0.6, sym='', widths=0.4)
bpm = ax1.boxplot(Data2, positions=np.array(range(len(Data2)))*3.0, sym='', widths=0.4)
bpr = ax1.boxplot(Data3, positions=np.array(range(len(Data3)))*3.0+0.6, sym='', widths=0.4)
set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
set_box_color(bpm, '#2C7BB6')
set_box_color(bpr, '#5D7FF6')

# draw temporary red and blue lines and use them to create a legend
ax1.plot([], c='#D7191C', label='Marine')
ax1.plot([], c='#2C7BB6', label='Original')
ax1.plot([], c='#5D7FF6', label='Introduced')

#ax1.set_title('Mean FAA per Individual Throughout Simulation')
ax1.set_ylabel('Mean number of alleles / individual',fontsize=10)
ax1.set_xlabel('$M$',fontsize=10)
ax1.set_title("(A)")
plt.xticks(range(0, len(ticks) * 3, 3), mig_rates)
plt.xlim(-1.0, len(ticks)*2.7)
#plt.tight_layout()

bpl = ax2.boxplot(DataP1, positions=np.array(range(len(DataP1)))*3.0-0.6, sym='', widths=0.4)
bpm = ax2.boxplot(DataP2, positions=np.array(range(len(DataP2)))*3.0, sym='', widths=0.4)
bpr = ax2.boxplot(DataP3, positions=np.array(range(len(DataP3)))*3.0+0.6, sym='', widths=0.4)
set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
set_box_color(bpm, '#2C7BB6')
set_box_color(bpr, '#5D7FF6')

ax2.set_ylabel('Mean percentage of alleles / ind',fontsize=10)
ax2.set_xlabel('$M$',fontsize=10)
ax2.set_title("(B)")

plt.xticks(range(0, len(ticks) * 3, 3), mig_rates)
plt.xlim(-1.0, len(ticks)*2.7)

bpl = ax3.boxplot(DataT, positions=np.array(range(len(DataT)))*3.0, sym='', widths=0.4)
set_box_color(bpl, '#000000')

ax3.set_ylabel('Total number of alleles defined',fontsize=10)
ax3.set_xlabel('$M$',fontsize=10)
ax3.set_title("(C)")

mpatch = mpatches.Patch(color='#D7191C',label='Marine')
cpatch = mpatches.Patch(color='#2C7BB6',label='Old lakes')
ypatch = mpatches.Patch(color='#5D7FF6',label='New Lakes')
ax1.legend(handles=[mpatch,cpatch,ypatch],loc = 1,frameon=False,fontsize=7)

plt.xticks(range(0, len(ticks) * 3, 3), mig_rates,label = '$m$ (migrants / lake / generation)')
plt.xlim(-1.0, len(ticks)*2.7)
fig.subplots_adjust(left=.10, bottom=.15, right=.98, top=.90,wspace = 0.4)
height = 6.00
width = 3.25
fig.set_size_inches(height, width)
fig.savefig("Freshwater_Alleles.pdf")

