import pandas as pd
import numpy as np
import glob
import os

path =r'/gpfs/home/goliveir/2021.01.22'
allFiles = glob.glob(path + "/*.tsv")
genomepos = np.arange(1, 10500, 1)
#cols = ['ALT_FREQ']
varfreqdf = pd.DataFrame(index=genomepos)
df = pd.DataFrame(index=genomepos)
dflist = []
files = {}
dflist.append(varfreqdf)

for varfile in allFiles:
	print(varfile)
	file_contents = open(varfile, 'r').read()
	files[varfile] = file_contents
	df = pd.read_csv(varfile,index_col=1, header=0,sep='\t')
	s = df.iloc[:,9]
	s= s.loc[~s.index.duplicated(keep='first')]
	dflist.append(s)
	
allFiles = [e[28:-8] for e in allFiles]
varfreqdf = pd.concat(dflist, axis=1)
print(varfreqdf)
varfreqdf.columns = allFiles
result = varfreqdf.iloc[[3069,3533]] 
result.to_csv('ZikaSingleAmpplate1.csv',sep='\t', encoding='utf-8')


print(allFiles)



