
import matplotlib as mpl
mpl.use('pdf')
import msprime
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

mig_rate = ["$0.01$","$0.1$","$1.0$","$10.0$","$100.0$"]
combinations = ["$F_{st}$"+"\nMarine\nOld Lakes","$F_{st}$"+"\nMarine\nNew Lakes","$F_{st}$"+"\nOld Lakes\nNew Lakes"]
x = [4949999, 5049999, 9999999, 14949999, 15049999, 19999999, 24949999, 25049999, 29999999, 34949999, 35049999, 39999999, 44949999, 45049999, 49999999, 54949999, 55049999, 59999999, 64949999, 65049999, 69999999, 74949999, 75049999, 79999999, 84949999, 85049999, 89999999, 94949999, 95049999, 99999999]

effectIntervals = []

for i in range(len(x)-1):
    if(i % 3 == 0):
        effectIntervals.append([x[i],x[i+1]])    

meanEffect = []
for e in effectIntervals:
    meanEffect.append(sum(e)/len(e))

treeSequences = []

def populateTreeSequences():

    for i in range(5):
        ts = msprime.load("./Output1/MyRecipe9_1_"+str(i)+"/EndTrees.trees") 
        mut_ts = msprime.mutate(ts,rate=1e-8,keep=True,random_seed=25)
        treeSequences.append(mut_ts)


def mutations_fst(tree_sequence, pop0, pop1, windows):
    '''
    Compute mean Fst across all variants in each window.
    Assumes an infinite-sites model.
    '''
    n0 = len(pop0)
    n1 = len(pop1)
    num_windows = len(windows) - 1
    if windows[0] != 0.0:
        raise ValueError(
            "Windows must start at the start of the sequence (at 0.0).")
    if windows[-1] != tree_sequence.sequence_length:
        raise ValueError("Windows must extend to the end of the sequence.")
    for k in range(num_windows):
        if windows[k + 1] <= windows[k]:
            raise ValueError("Windows must be increasing.")
    S = np.zeros(num_windows)
    # index of *left-hand* end of the current window
    window_num = 0
    num_sites = 0
    tt0 = tree_sequence.trees(tracked_samples=pop0,
                                   sample_counts=True)
    tt1 = tree_sequence.trees(tracked_samples=pop1,
                                   sample_counts=True)
    for t0, t1 in zip(tt0, tt1):
        for s in t0.sites():
            while s.position >= windows[window_num + 1]:
                if num_sites > 0:
                    S[window_num] /= num_sites
                window_num += 1
                num_sites = 0
            num_sites += 1
            for m in s.mutations:
                x0 = t0.num_tracked_samples(m.node)
                x1 = t1.num_tracked_samples(m.node)
                p0 = x0 / n0
                p1 = x1 / n1
                if (x0 + x1 > 0) and (x0 + x1 < n0 + n1):
                    fst = 1 - (p0 * (1 - p0) + p1 * (1 - p1)) / ((p0 + p1) * (1 - (p0 + p1) / 2))
                    S[window_num] += fst
    # do the last window
    if num_sites > 0:
        S[window_num] /= num_sites
    return S

def inEffectRegion(windows,windowSize):

    inEffect = np.zeros(len(windows)-1)
    
    for k in range(len(effectIntervals)):
        a = effectIntervals[k][0]
        b = effectIntervals[k][1]
        leftmostWindow = a//windowSize
        rightmostWindow = b//windowSize
        for win in range(leftmostWindow,rightmostWindow):
            inEffect[win] = 1
        inEffect[rightmostWindow] = 1
        
    return inEffect

def showSpatialLocations(tree_sequence,pop):

    inds = tree_sequence.individuals()

    y = []
    x = []
    count = 0 
    
    for i in inds:
        node = ts.node(i.nodes[0])
        if(node.population == pop):
            y.append(i.location[0])
            x.append(count)
            count += 1
    plt.scatter(x,y)
    plt.show()    

def getFst(tree_sequence,left_location_bound,right_location_bound,pop0,pop1,window_size):

    windows = [i for i in range(0,int(tree_sequence.sequence_length),window_size)] + [int(tree_sequence.sequence_length)]
    inds = tree_sequence.individuals()
    p1 = []
    p2 = []

    for i in inds:
        if((i.location[0] >= left_location_bound) and (i.location[0] < right_location_bound)):
            node0_id = i.nodes[0]
            node1_id = i.nodes[1]
            node0 = tree_sequence.node(node0_id)
            node1 = tree_sequence.node(node1_id)
            if(node0.population == pop0):
                p1.append(node0_id)
                p1.append(node1_id)
            elif(node0.population == pop1):
                p2.append(node0_id)
                p2.append(node1_id)

    fst1_2 = mutations_fst(tree_sequence = tree_sequence, pop0 = p1, pop1 = p2, windows = windows) 
    windows = np.array(windows)
    fst1_2 = np.array(fst1_2)
    return windows,fst1_2
    
def WriteFsts():

    for k,combs in enumerate([[1,2],[1,3],[2,3]]):
        for i in range(len(mig_rate[:-1])):

            mut_ts = treeSequences[i]

            x,y = getFst(tree_sequence = mut_ts,
                        left_location_bound = 0.0,
                        right_location_bound = 25.0,
                        pop0 = combs[0],
                        pop1 = combs[1],
                        window_size = 250
                        )
            x = " ".join(x.astype(str))            
            y = " ".join(y.astype(str))
            strings = [x,"\n",y]    

            path = "./fsts_250"
            filename = os.path.join(path,str(k)+"_"+str(i)+"_fst.txt")
            fi = open(filename,"w")
            fi.writelines(strings)

    return None


opt = "0_5_500"
#opt = ""

def ReadAndPlotFsts():

    fig, axes = plt.subplots(nrows=3, ncols=len(mig_rate[:-1]), figsize=(6, 6), sharex=True, sharey=True)
    colors = ["c","k","y"]
    
    prefix_faa_loc = "FAA_Positions/9_1_"
    postfix_faa_loc = "FAA_Positions.txt"
    handles = []

    for k in range(3):
        for i in range(4):
            path = "./fsts" + opt
            filename = os.path.join(path,str(k)+"_"+str(i)+"_fst.txt")
            fi = open(filename)
            x = [float(X) for X in fi.readline().split()]
            y = [float(X) for X in fi.readline().split()]

            faa_loc_file = open(prefix_faa_loc+str(i)+postfix_faa_loc,"r")
            positions = faa_loc_file.readline().split()
            for p in positions:
                axes[k][i].axvline(int(p),linestyle = '--',linewidth=0.3, color='pink',alpha=0.2)
            line = axes[k,i].plot(x[1:], y, lw = 0.45, c=colors[k])
            if(k == 0):
                axes[k][i].set_title(mig_rate[i])
            if(i == 0):
                axes[k][i].set_ylabel(combinations[k])
    
    fig.subplots_adjust(left=.18, bottom=.15, right=.95, top=.90)
    height = 3.75
    width = 6.00
    fig.set_size_inches(width, height)
    fig.savefig("Fst_Genome_faa"+opt+".pdf")

def TrueNegative():

    TNs = []
    SPs = []
    FSTs = []

    for i in range(len(mig_rate)-1):
        
        path = "./fsts_250"
        filename = os.path.join(path,"1_"+str(i)+"_fst.txt")
        fi = open(filename)
        x = np.array([float(X) for X in fi.readline().split()])
        y = np.array([float(X) for X in fi.readline().split()])
        
        inEffect = inEffectRegion(x,250)

        TN = np.cumsum(inEffect[y.argsort()[::-1]]) / np.arange(1,len(inEffect)+1)
        SP = np.cumsum(inEffect[y.argsort()[::-1]]) / inEffect.sum()
        y.sort()
        TNs.append(TN)
        SPs.append(SP)
        FSTs.append(y[::-1])

    return TNs,SPs,FSTs

def plotTrueNegative():

    TNs,SPs,FSTs = TrueNegative()        

    fig1,ax1 = plt.subplots(2,1,sharex = True)

    for k in range(len(mig_rate)-1):
        ax1[0].plot(FSTs[k],TNs[k],label = mig_rate[k])
    for k in range(len(mig_rate)-1):
        ax1[1].plot(FSTs[k],SPs[k],label = mig_rate[k])

    ax1[0].set_title('True Positive')
    ax1[0].grid()
    ax1[0].legend(loc='upper right',frameon=False)
    ax1[1].set_title('Statistical Power')
    ax1[1].grid()

    fig1.subplots_adjust(left=.15, bottom=.10, right=.85, top=.90)
    height = 7.00
    width = 3.75
    fig1.set_size_inches(height, width)
    fig1.savefig("True_Power_250.pdf")


if __name__ == "__main__":
    

    # For plotting Fst across genome, First run
    # populateTreeSequences()
    # WriteFsts()
    # Then after directories are created, run
    # ReadAndPlotFsts()

    # For plotting StatPower, run
    # populateTreeSequences()
    # plotTrueNegative()
    

    




