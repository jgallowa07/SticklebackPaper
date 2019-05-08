import matplotlib as mpl
mpl.use('pdf')
import msprime
import numpy as np
import pyslim
import matplotlib.pyplot as plt

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

samples_prop_de_novo = []
samples_prop_migrant = []
samples_prop_standing = []


def mutationsLocation(mutationSlimIds,treeSequence):
    '''
    take in SLiM mutation Id's and return their respective location using pyslim.
    '''

    positions = np.empty(len(mutationSlimIds),dtype=int)
    mutation_msp_ids = np.empty(len(mutationSlimIds),dtype=int)
    mutations = treeSequence.mutations()
    found_count = 0
    for m in mutations:
        try:
            idx = mutationSlimIds.index(int(m.derived_state))
            positions[idx] = int(m.position)
            mutation_msp_ids[idx] = int(m.id)
            found_count += 1
        except:
            continue
    #print(found_count," ",len(mutationSlimIds))
    #print(positions)

    assert(found_count == len(mutationSlimIds))
    return positions, mutation_msp_ids

def countFitnessOfAncestors(mutationSlimIds,treeSequence):

    positions,mutationMspIds = mutationsLocation(mutationSlimIds,treeSequence)

    #Get the extant samples from p3 and the initial population from p3. 
    #They are all marked as samples in the treeSequence
    extantSampleIds = []
    IntroducedPopulationsIds = [] 
    for sam in treeSequence.samples():
        node = treeSequence.node(sam)
        time = node.time
        pop = node.population
        if time == 0.0 and pop == 3:
            extantSampleIds.append(sam)
        elif pop == 3:
            IntroducedPopulationsIds.append(sam)
            T = time
        else:
            continue

    fitnesses = np.array([0 for _ in IntroducedPopulationsIds])
    #print(len(extantSampleIds))
    #return
    tt = ts.trees(tracked_samples=extantSampleIds)
    t = next(tt)
    for var in ts.variants():
        if (var.site.position in positions):
            while t.interval[0] <= var.site.position:
                t = next(tt)
            for i,k in enumerate(IntroducedPopulationsIds):
                fitnesses[i] += t.num_tracked_samples(k)
        else:
            continue
    return fitnesses    

def originOfAdaptiveVariants_dev(mutationSlimIds,preExistingMutationSlimIds,treeSequence):

    positions,mutationMspIds = mutationsLocation(mutationSlimIds,treeSequence)

    extantSampleIds = []
    IntroducedPopulationsIds = [] 
    for sam in treeSequence.samples():
        node = treeSequence.node(sam)
        time = node.time
        pop = node.population
        if time == 0.0 and pop == 3:
            extantSampleIds.append(sam)
        elif pop == 3:
            IntroducedPopulationsIds.append(sam)
            T = time
        else:
            continue

    treeStats = []
    originNodeCounts = []

    tt = ts.trees(tracked_samples=extantSampleIds)
    t = next(tt)
    for var in ts.variants():
        if (var.site.position in positions):
            while t.interval[0] <= var.site.position:
                t = next(tt)
        else:
            continue
    return fitnesses    

def originOfAdaptiveVariants(mutationSlimIds,preExistingMutationSlimIds,treeSequence):
    '''
    build the trees at locations of adaptive variants and
    ask, how many extant samples have the introduced populaton as MRCA at this site. 

    this will return a tuple for each tree that a freshwater adapted allele falls on
    the indices will contain a count of the number of samples who have inherited that region of the genome from a
    migrant, original introduced population, or have a deNovo mutation, respectively.  
   
    '''
    positions,mutationMspIds = mutationsLocation(mutationSlimIds,treeSequence)

    prePositions,preMutationMspIds = mutationsLocation(preExistingMutationSlimIds,treeSequence)

    #Get the extant samples from p3 and the initial population from p3. 
    #They are all marked as samples in the treeSequence
    extantSampleIds = []
    IntroducedPopulationsIds = [] 
    for sam in treeSequence.samples():
        node = treeSequence.node(sam)
        time = node.time
        pop = node.population
        if time == 0.0 and pop == 3:
            extantSampleIds.append(sam)
        elif pop == 3:
            IntroducedPopulationsIds.append(sam)
            T = time
        else:
            continue


    treeStats = []
    originNodeCounts = []
    
    for tree in treeSequence.trees():
        
        interval = tree.interval
        positionInQuestion = 0 
        mutationInQuestion = 0
        mutationInQuestionIsStanding = False
        pos_count = 0
        for i in range(len(positions)):
            if(positions[i] >= interval[0] and positions[i] < interval[1]):
                positionInQuestion = positions[i]
                mutationInQuestion = mutationMspIds[i] 
                if(mutationInQuestion in preMutationMspIds):
                    mutationInQuestionIsStanding = True
                pos_count += 1
        if (pos_count >= 1):
            samples,deNovo = samplesWithAdaptiveVariant(treeSequence,tree,extantSampleIds,mutationInQuestion)
            origin_migrants = 0
            origin_standing = 0
            origin_captured = 0  
            origin_deNovo = 0
            standingIds = []

            #If the mutation in question is deNovo, then all samples who have the 
            #mutation should count towards the deNovo category.
            if(deNovo):
                treeStats.append(tuple([origin_migrants,origin_standing,origin_captured,len(samples)]))    
                originNodeCounts.append(0)
            else:
                for s in samples:
                    Nod = lastLake(s,tree,treeSequence,T)
                    if(Nod in IntroducedPopulationsIds):
                        if(mutationInQuestionIsStanding):
                            origin_captured += 1
                        else:
                            origin_standing += 1
                        standingIds.append(Nod)
                    else:
                        origin_migrants += 1

                treeStats.append(tuple([origin_migrants,origin_standing,origin_captured,origin_deNovo]))
                originNodeCounts.append(len(set(standingIds)))

    return treeStats,originNodeCounts
                

def lastLake(u,sparseTree,treeSequence,T):

    v = sparseTree.parent(u)
    while (v != msprime.NULL_NODE and treeSequence.node(v).population == 3 and treeSequence.node(v).time <= T):
        u = v
        v = sparseTree.parent(u)
    return u 

def lastLakeAverageBL(u,sparseTree,treeSequence,T):

    BranchLengths = []
    v = sparseTree.parent(u)
    while (v != msprime.NULL_NODE and treeSequence.node(v).population == 3 and treeSequence.node(v).time <= T):
        BranchLengths.append(sparseTree.branch_length(u))
        u = v
        v = sparseTree.parent(u)
    averageBL = sum(BranchLengths)/len(BranchLengths)    

    return averageBL

def samplesWithAdaptiveVariant(TreeSequence,sparseTree,samples,mutationId):
    deNovo = False
    mutation = TreeSequence.mutation(mutationId)
    originNodeOfMutation = mutation.node
    pop = TreeSequence.node(originNodeOfMutation).population
    if(pop == 3):
        deNovo = True
    samplesWithMutation = []
    for s in samples:
        mrca = sparseTree.mrca(s,originNodeOfMutation)
        if mrca == originNodeOfMutation:
            samplesWithMutation.append(s)
        
    return samplesWithMutation,deNovo    
 

def readFiles(recipeNum):
    '''
    read in the tree sequence and the mutations Ids to find origin
    '''


    adaptationFile = open("./Output1/MyRecipe"+recipeNum+"/Adaptation.txt","r")
    ts = pyslim.load("./Output1/MyRecipe"+recipeNum+"/TreesAtAdaptation.trees",slim_format=True)

    header = adaptationFile.readline()
    PreExisting_FWAA = [int(Id) for Id in adaptationFile.readline().split()]
    Introduced_FWAA = [int(Id) for Id in adaptationFile.readline().split()]
    
    return Introduced_FWAA,PreExisting_FWAA,ts


'''
#COMPUTE AND PLOT FITNESSES
fits = []
numVariants = []
for i in range(4):
    fa,fa2,ts = readFiles("9_1_"+str(i))
    fit = countFitnessOfAncestors(fa,ts)
    numVariants.append(len(fa))
    fits.append(fit)

bins = [x for x in range(0,270,5)]

fig, axs = plt.subplots(4,sharex=True)

for i,k in enumerate(fits):
    axs[i].hist(k[k!=0]/numVariants[i],bins=bins)

plt.show()

#COMPUTE AND PLOT FITNESSES
fits = []
numVariants = []
for i in range(4):
    fa,fa2,ts = readFiles("9_1_"+str(i))
    fit = averageBranchLengthsOfVariants(fa,ts)
    #numVariants.append(len(fa))
    fits.append(fit)

#bins = [x for x in range(0,270,5)]

fig, axs = plt.subplots(4,sharex=True)

for i,k in enumerate(fits):
    axs[i].hist(k[k!=0])

plt.show()
'''
'''
treeStats = []
sample_counts = []
for i in range(4):
    fa_p3,fa_p2,ts = readFiles("9_1_"+str(i))
    treeStat,counts = originOfAdaptiveVariants(fa_p3,fa_p2,ts)
    treeStats.append(treeStat)
    sample_counts.append(counts)


mig = []
pre = []
cap = []
deN = []

for mig_rate in treeStats:
    mig.append(sum([i[0] for i in mig_rate]))
    pre.append(sum([i[1] for i in mig_rate]))
    cap.append(sum([i[2] for i in mig_rate]))
    deN.append(sum([i[3] for i in mig_rate]))

print(mig)
print(pre)
print(cap)
print(deN)
'''

#mig = [2379, 2847, 1314, 5881]
#pre = [20656, 83912, 62169, 44220]
#deN = [89633, 3357, 0, 0]

mig = [2379, 2847, 1314, 5881]
pre = [1656, 3000, 0, 0]
cap = [19000, 80912, 62169, 44220]
deN = [89633, 3357, 0, 0]

ratio_mig = []
ratio_pre = []
ratio_cap = []
ratio_deN = []

for i in range(len(mig)):
    s = mig[i] + pre[i] + cap[i] + deN[i] 
    ratio_mig.append(mig[i] / s)
    ratio_pre.append(pre[i] / s)
    ratio_cap.append(cap[i] / s)
    ratio_deN.append(deN[i] / s)

ind = np.arange(4)
width = 0.35

cap_bottom = []
for i in range(len(mig)):
    cap_bottom.append(mig[i] + pre[i])

den_bottom = []
for i in range(len(mig)):
    den_bottom.append(cap_bottom[i] + cap[i])

ratio_cap_bottom = []
for i in range(len(ratio_mig)):
    ratio_cap_bottom.append(ratio_mig[i] + ratio_pre[i])

ratio_den_bottom = []
for i in range(len(mig)):
    ratio_den_bottom.append(ratio_cap_bottom[i] + ratio_cap[i])

fig,axes = plt.subplots()

p1 = axes.bar(ind, ratio_mig, width, color = "#5780CD")
p2 = axes.bar(ind, ratio_pre, width, bottom=ratio_mig,color="#F06261")
p3 = axes.bar(ind, ratio_cap, width, bottom=ratio_cap_bottom,color="#A958A5")
p4 = axes.bar(ind, ratio_deN, width, bottom=ratio_den_bottom,color="#42B97C")

plt.ylabel('Allele origin ratio')
plt.xlabel('$M$')
plt.xticks(ind, ('0.01', '0.1', '1.0', '10.0'))
plt.legend((p1[0], p2[0], p3[0], p4[0]), 
            ('Migrants','Marine','Captured','De novo'),
            loc = 'upper center',bbox_to_anchor=(1.15,1.0),frameon=False
            ,fontsize=10
            )

fig.subplots_adjust(left=.15, bottom=.16, right=.80, top=.97)
width = 6.0
height =  3.0
fig.set_size_inches(width, height)
fig.savefig("./figs/Allele_Origin_2.pdf")
