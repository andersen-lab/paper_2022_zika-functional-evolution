
import glob
import textwrap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg') 	
import matplotlib.pyplot as plt
from scipy.special import xlogy
import seaborn as sns
import scipy.stats as stats
aminocode = {"ATG":'Start',"GCT":'A',"GCC":'A',"GCA":'A',"GCG":'A',"GGT":'G',"GGC":'G',"GGA":'G',"GGG":'G',"ATT":'I',"ATC":'I',"ATA":'I',"CTT":'L',"CTC":'L',"CTA":'L',"CTG":'L',"TTA":'L',"TTG":'L',"CCT":'P',"CCC":'P',"CCA":'P',"CCG":'P',"GTT":'V',"GTC":'V',"GTA":'V',"GTG":'V',"TTT":'F',"TTC":'F',"TGG":'W',"TAT":'Y',"TAC":'Y',"GAT":'D',"GAC":'D',"GAA":'E',"GAG":'E',"CGT":'R',"CGC":'R',"CGA":'R',"CGG":'R',"AGA":'R',"AGG":'R',"CAT":'H',"CAC":'H',"AAA":'K',"AAG":'K',"TCT":'S',"TCC":'S',"TCA":'S',"TCG":'S',"AGT":'S',"AGC":'S',"ACT":'T',"ACC":'T',"ACA":'T',"ACG":'T',"TGT":'C',"TGC":'C',"ATG":'M',"AAT":'N',"AAC":'N',"CAA":'Q',"CAG":'Q',"TGA":'X',"TAA":'X',"TAG":'X'}
clades = [2074 , 3045,894 , 2787,1143 , 2842,3398 , 3328,107 , 1118,3353 , 80,139 , 2086,2634 , 988, 1542]
content = []
split = []
muts = []
muts2 = []
mutsdnds = []
df = pd.DataFrame()
mutdf = pd.DataFrame(columns=["Posit", "RefNT","Mutname", "Mut","MutCodon","RefCodon","Partition","RefAA","MutAA","Type"])
with open('TimeTreeAncRestNode1_012420.fa', 'r') as file:
	Ref = file.read().replace('\n', '')
ZikaRef = Ref[13:]
nt_list = [c for c in ZikaRef]

df['RefNT'] = nt_list
df.index = np.arange( 1, len(df) + 1)


AAlist = textwrap.wrap(ZikaRef, 3)

AAlist = list(np.repeat(AAlist, 3))

print(AAlist)
print(df)
print(len(AAlist))
print(len(df))
df['RefCodon'] = AAlist


df['RefAA'] = df['RefCodon'].map(aminocode)
#pd.set_option('display.max_rows', None)



import re
with open('/Users/glennoliveira/Documents/2021-06-24-0003_ancestral/annotated_tree.nexus','r') as f:
    for ln in f:
    	content.append(re.findall(r'"(.*?)"', ln))
flat_list = [item for sublist in content for item in sublist]

for s in flat_list:
	s = s.split(',')
	print(s)
	muts.extend(s)

#remove deletions
muts = [x for x in muts if not '-' in x]
print(muts)
ntlist = ['A','T','C','G','a','t','c','g']
countix = 0 
for i in muts:
	if i != "" and i[-1] in ntlist:
		mutname = i
		print(i)
		#mut.append(i)
		#Get rid of AA on both sides of locus
		mutnt = i[-1]
		print(mutnt)
		
		i = i[1:len(i)-1]
		i = int(i)# -1
		print(i)
		mutdf = mutdf.append(df.loc[i],ignore_index=True)
		mutdf.loc[countix,'Mut'] = mutnt
		mutdf.loc[countix,'Mutname'] = mutname
		mutdf.loc[countix,'Posit'] = i
	
		if int(i) % 3 == 0:
			#print('multof3')
			i = i/3
			#print(i)
			mutdf.loc[countix,'Partition'] = 2
		elif int(i) % 3 == 1:
			i = (i-1)/3
			#print('1left')
			#print(i)
			mutdf.loc[countix,'Partition'] = 0
		elif int(i) % 3 == 2:
			i = (i-2)/3
			#print('2left')
			#print(i)
			mutdf.loc[countix,'Partition'] = 1
		AAseq = list(str(mutdf.loc[countix,'RefCodon']))
		mutnt = str(mutdf.loc[countix,'Mut'])
		AApos = int(mutdf.loc[countix,'Partition'])
		
	#AAseq[mutdf.loc[countix,'Partition']] = mutdf.loc[countix,'Mut']
	
	AAseq[AApos] = mutnt

	
	mutAAseq = "".join(AAseq)
	
	mutdf.loc[countix,'MutCodon'] = mutAAseq
	

	countix = countix +1 
mutdf['MutAA'] = mutdf['MutCodon'].map(aminocode)

mutdf["Type"] = np.where(mutdf['RefAA'] == mutdf['MutAA'],'S', 'NS')
mutdf['AAposit'] = (mutdf['Posit'] + 2) // 3
                                           
print(mutdf.Type.value_counts())
print(mutdf.sort_values(by=['AAposit']))
print(mutdf)

df = pd.read_csv('logcombined062221dNdSFinalsampledat35k.log',skiprows=4,index_col=0, sep='\t',header=None)
df.loc['mean'] = df.mean()
mean = df.loc['mean']
#clades = mean.ix[[2074 , 3045,894 , 2787,1143 , 2842,3398 , 3328,107 , 1118,3353 , 80,139 , 2086,2634 , 988]]
#mean = mean.drop([2074 , 3045,894 , 2787,1143 , 2842,3398 , 3328,107 , 1118,3353 , 80,139 , 2086,2634 , 988])
mean = mean.to_frame()
mean['AAposit'] = mean.index

finaldf = mutdf.merge(mean,left_on=['AAposit'], right_on=['AAposit'], how='left')
print(finaldf.loc[finaldf["Type"] == "NS"])
finaldf = finaldf.drop_duplicates(subset=['AAposit'])




finaldf.loc[finaldf["Type"] == "NS"].to_csv('ZFEnonsyn.06.29.21.csv')
finaldf.loc[finaldf["Type"] == "S"].to_csv('ZFEsynon.06.29.21.csv')

#print(finaldf)
#print(finaldf.sort_values(by=['AAposit']))
for c in clades:
	mask = finaldf['AAposit'] == c
	finaldf['Type'][mask] = 'C' 

finaldf.loc[finaldf["Type"] == "C"].to_csv('ZFEclades.06.29.21.csv')

#print(finaldf.sort_values(by=['mean']))
#print(finaldf)
clades = finaldf.loc[finaldf["Type"] == "C"]
#print(clades)
#print(finaldf.loc[finaldf["Type"] == "NS"])
clades = clades.drop_duplicates(subset=['AAposit'])

nonsynmuts = finaldf.loc[finaldf["Type"] == "NS"]
synmuts = finaldf.loc[finaldf["Type"] == "S"]
#nonsynmuts = nonsynmuts.drop_duplicates(subset=['AAposit'])
u_statistic, pVal = stats.mannwhitneyu(clades['mean'], nonsynmuts['mean'])

print(nonsynmuts)
print(clades)

print("Here is the pval:")

print(pVal)

boxdf = pd.DataFrame()

synmuts = synmuts['mean'].dropna()
nonsynmuts = nonsynmuts['mean'].dropna()
clades = clades['mean'].dropna()


boxdf = pd.concat([synmuts,nonsynmuts,clades], axis=1,)
boxdf.columns = ['Lineage-defining','Nonsynonymous','Synonymous']



u_statistic, pVal = stats.mannwhitneyu(clades, nonsynmuts)

print(nonsynmuts)
print(clades)
print("Here is the pval:")

print(pVal)
font = {'family' : 'normal',
        'size'   : 30}

matplotlib.rc('font', **font)

plt.figure(figsize=(17,8))
ax = plt.boxplot([synmuts,nonsynmuts,clades], labels=['Synonymous','Nonsynonymous','Lineage-defining'],showfliers=False,widths=.75)



for i, d in enumerate(boxdf):
   y = boxdf[d]
   x = np.random.normal(i + 1, 0.1, len(y))
   plt.scatter(x, y,s=190,c='royalblue',edgecolors='black')
plt.axhline(y=1, color='black', linestyle='dashed', linewidth=1)

plt.savefig('ZFEmutsdndstatsboxandwhisker2.pdf', dpi=300,format='pdf')

