import matplotlib as mpl
mpl.use('pdf')
import msprime
import pyslim
import io
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

ts = pyslim.load("./Output1/MyRecipe9_1_2/EndTrees.trees")

inds = [ind.id for ind in ts.individuals()]
sample_subset = np.unique(np.random.choice(inds,6000))
sample_nodes = []
for i in sample_subset:
    ind = ts.individual(i)
    sample_nodes.append(ind.nodes[0])
    sample_nodes.append(ind.nodes[1])
    
ts = ts.simplify(sample_nodes,filter_populations=False)

extant_samples = []
for n in ts.nodes():
    if (n.time == 0.0):
        extant_samples.append(n.id)

tables = ts.dump_tables()
x = [[],[],[]]
y = [[],[],[]]
c = [[],[],[]]
varit = ts.variants()
colors = ['m','c','y']

ind_loc = np.array([i.location[0] + np.random.uniform(-0.5,0.5) for i in ts.individuals()])
phenotypes = np.array([0.0 for _ in ts.individuals()])
populations = np.array([ts.node(ind.nodes[0]).population for ind in ts.individuals()])

for v in varit:
    mut_meta = pyslim.decode_mutation(v.site.mutations[0].metadata)[0]
    selectionCoeff = mut_meta.selection_coeff
    mutType = mut_meta.mutation_type
    color = None
    dominance = None
    if(mutType == 1 or mutType == 2):
        color = 'm'
        dominance = 'dom'
    elif(mutType == 4 or mutType == 5):
        color = 'c'
        dominance = 'rec'
    else:
        color = 'y'
        dominance = 'add'

    for ind in ts.individuals():
        num_muts = sum([v.genotypes[n] for n in ind.nodes])
        if dominance == 'dom' and num_muts > 0:
            phenotypes[ind.id] += 2 * selectionCoeff;
        if dominance == 'rec' and num_muts == 1:
            phenotypes[ind.id] += 0;
        if dominance == 'rec' and num_muts == 2:
            phenotypes[ind.id] += 2 * selectionCoeff;
        if dominance == 'add':
            phenotypes[ind.id] += num_muts * selectionCoeff;
    
    for i,g in enumerate(v.genotypes):
        if (g == 1 and i in extant_samples):
            node = ts.node(i)
            pop = node.population
            sp = ind_loc[node.individual]
            x[pop-1].append(sp)
            y[pop-1].append(selectionCoeff)
            c[pop-1].append(color)

fig,axes = plt.subplots(nrows = 2,ncols = 3,sharex = True,sharey = 'row')
labels = ['Marine','Old Lakes','New Lakes']
for l in range(3):
    axes[0,l].scatter(x=ind_loc[populations==(l+1)],y=phenotypes[populations==(l+1)],marker='.',s=0.1)
    axes[0,l].set_title(labels[l])
    if (l == 0):
        axes[0,l].set_ylabel('Trait value')
    axes[0,l].set_ylim(-15,13)
for k in range(3):
    axes[1,k].scatter(x=x[k],y=y[k],c=c[k],marker='.',s=0.1,alpha=0.7,label=labels[l])
    if(k == 0):
        axes[1,k].set_ylabel('Effect size')
    axes[1,k].set_xlabel('Spatial position')

mpatch = mpatches.Patch(color='m',label='dominant')
cpatch = mpatches.Patch(color='c',label='recessive')
ypatch = mpatches.Patch(color='y',label='additive')

plt.legend(handles=[mpatch,cpatch,ypatch],loc='center left', bbox_to_anchor=(1, 0.5),frameon = False)
fig.subplots_adjust(left=.15, bottom=.15, right=.80, top=.90)
width = 6.00
height = 3.75
fig.set_size_inches(width, height)
fig.savefig("Haplo_small.pdf")

