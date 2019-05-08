import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import sys
import numpy as np

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

def readHeader(firstLine):

	numberOfSamplingPoints,samplingInterval,introduction,endofSim = (int(i) for i in firstLine.split())
	samplesUntilIntroduction = introduction//samplingInterval - 1
	return numberOfSamplingPoints,samplingInterval,introduction,endofSim,samplesUntilIntroduction

def phenoData(recipeNum):

    data = []
    numLakes = 25

    for i in range(0,((2*numLakes)+3)):
        data.append([])

    x_axis = []
    x_axis2 = []
    lake_traces = []
    File = open("./Output1/MyRecipe"+recipeNum +"/AveragePhenotypeThroughout.txt","r")
    first = File.readline()
    numberOfSamplingPoints,samplingInterval,introduction,end,samplesUntilIntroduction = readHeader(first)
    print(numberOfSamplingPoints,samplingInterval,introduction,end,samplesUntilIntroduction)

    for i,line in enumerate(File):
        for j,s in enumerate(line.split()):
            data[j].append(float(s))

    #Give all data a zero point because that wasn't sampled
    for i in range(0,(numLakes+2)):
        data[i] = [0] + data[i]

    for i in range(0,numberOfSamplingPoints + 1):
        x_axis.append(i*samplingInterval)

    for i in range(samplesUntilIntroduction,numberOfSamplingPoints):
        x_axis2.append(i*samplingInterval)
    
    return data, x_axis, x_axis2

firstdata,x1,x2 = phenoData("9_1_0")
seconddata,x11,x22 = phenoData("9_1_1")

fig,axes = plt.subplots(nrows=2,ncols=1,sharex = True,sharey = True)

lw1 = 0.7
lw2 = 0.03


axes[0].plot(x1,firstdata[0],label='Marine',lw=lw1)
axes[0].plot(x1,firstdata[1],label='Old lakes',lw=lw1)
axes[0].plot(x2,firstdata[27],label='New lakes',lw=lw1)

for i in range(0,25):    
    axes[0].plot(x1,firstdata[i+2],c="y",lw=lw2)

for i in range(0,25):
    axes[0].plot(x2,firstdata[i+27],c='c',lw=lw2)

axes[1].plot(x11,seconddata[0],lw=lw1)
axes[1].plot(x11,seconddata[1],lw=lw1)
axes[1].plot(x22,seconddata[27],lw=lw1)

for i in range(0,25):    
    axes[1].plot(x11,seconddata[i+2],c="y",lw=lw2)

for i in range(0,25):
    axes[1].plot(x22,seconddata[i+27],c='c',lw=lw2)

axes[0].set_ylabel("Trait value")
axes[0].set_title("$(A)$")
axes[1].set_ylabel("Trait value")
axes[1].set_title("$(B)$")
axes[1].set_xlabel("Generations")

axes[0].legend(frameon=False,loc=0,fontsize = 8)


for i in range(2):
    for item in ([axes[i].xaxis.label, axes[i].yaxis.label] +
                 axes[i].get_xticklabels() + axes[i].get_yticklabels()):
        item.set_fontsize(12)

fig.subplots_adjust(left=.15, bottom=.16, right=.85, top=.92,hspace = 0.3)
height = 3.75
width = 6.00
fig.set_size_inches(width, height)
fig.savefig("./figs/Pheno_Time.pdf")





