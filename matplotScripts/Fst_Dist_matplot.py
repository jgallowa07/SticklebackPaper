import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys

num = sys.argv[1]
recipes = ['8_'+num+'_0','8_'+num+'_1','8_'+num+'_2','8_'+num+'_3']

def readHeader(firstLine):

	numberOfSamplingPoints,samplingInterval,introduction,endofSim = (int(i) for i in firstLine.split())
	samplesUntilIntroduction = introduction/samplingInterval - 1
	return numberOfSamplingPoints,samplingInterval,introduction

Original_Intro = []
Marine_Original = []
Marine_Intro = []

data = [[][][]]

	
for i in range(4):

	file1 = "../../Output1/MyRecipe"+recipes[i]+"/FreshwaterFreshwater2Fst.txt"
	#,"Freshwater","IntroducedFreshwater",'rgb(60,179,113)'
	file2 = "../../Output1/MyRecipe"+recipes[i]+"/OceanFreshwaterFst.txt"
	#,"Ocean","Freshwater",'rgb(127,255,212)'
	file3 = "../../Output1/MyRecipe"+recipes[i]+"/OceanFreshwater2Fst.txt"
	#,"Ocean","IntroducedFreshwater",'rgb(138,43,226)'

	

	for k,f in enumerate([file1,file2,file3]):

		File = open(f,"r")
		first = File.readline()

		#read in respective Fst/Pos pairs - they are in order.
		Fst = list(map(float,File.readline().split()))
		Pos = list(map(int,File.readline().split()))
		
		#sort both lists by Fst value
		Pos_Sorted = [fst for _,fst in sorted(zip(Fst,Pos))]
		Fst_Sorted = sorted(Fst)

		data[k].append()
		
	




