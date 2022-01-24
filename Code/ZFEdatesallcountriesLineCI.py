import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg') 	
from astropy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import datetime as dt
import matplotlib.dates as mdates
from matplotlib import gridspec
matplotlib.rc('ytick', labelsize=14) 
pd.set_option("display.max_rows", None, "display.max_columns", None)

df = pd.read_csv("2020.05.18.ZIKVfinal.tsv",sep="	",header=0,index_col=0)
print(df)
#df = df[df['Date'].str.len() >= 7]
print(df)
#df = df[df['country'] != 'PA']
#df['Month'] = df['Date'].str[:7]
df['Date'] = pd.to_datetime(df['Date'])
df['Date'] = df['Date'].dt.strftime('%Y-%m')


count_series = df.groupby(['Clade', 'Date']).size()
df = count_series.to_frame(name = 'size').reset_index()


colorsdict = {"CladeA":"#99CC33","CladeB":"#FFD43B", "CladeC":"#DAA520","CladeD":"#FF8C00","CladeE":"#3399FF","CladeF":"#9370DB", "CladeG":"#800080","CladeH":"#32CD32","CladeI":"#013220"}

d = {'Date': 'last' ,'CladeA': 'sum','CladeB': 'sum', 'CladeC': 'sum','CladeD': 'sum','CladeE': 'sum','CladeF': 'sum', 'CladeG': 'sum','CladeH': 'sum','CladeI': 'sum'}
fig, ax = plt.subplots(figsize=(10,7))  

clades = df['Clade'].drop_duplicates()

colors = ["#99CC33", "#FFD43B","#DAA520",'#FF8C00','#3399FF',"#9370DB",'#800080','#008000','#013220']
#df['Clade'].apply(lambda x: colors[x])
print(df)
df2 = pd.DataFrame(df['Date']).drop_duplicates()
df['Bottom'] = 0
pivot_df = df.pivot(index='Date', columns='Clade', values='size').fillna(0)

pivot_df1 = pivot_df.reset_index()


newdf = pivot_df.reset_index()
print(newdf)
df2015 = newdf.iloc[:18].sum(axis = 0)
df20161 = newdf.iloc[18:22].sum(axis = 0)
df20162 = newdf.iloc[22:26].sum(axis = 0)
df20163 = newdf.iloc[26:30].sum(axis = 0)
df2017 = newdf.iloc[30:].sum(axis = 0)
testdf = pd.concat([df2015,df20161,df20162,df20163,df2017], axis=1).T
testdf = testdf.drop(["Date",'PA1','PA2','PA3'], 1)


newdf = newdf.groupby(newdf.index // 11).agg(d)


newdf = newdf.set_index('Date')
newdf["sum"] = newdf.sum(axis=1)
testdf["sum"] = testdf.sum(axis=1)


trials = testdf["sum"] 
print(testdf)

cidf = testdf.iloc[: , :-1]

#get percentages
testdf[['CladeA','CladeB','CladeC','CladeD','CladeE','CladeF','CladeG','CladeH', 'CladeI']] = testdf[['CladeA','CladeB','CladeC','CladeD','CladeE','CladeF','CladeG','CladeH', 'CladeI']].div(testdf['sum'], axis=0)
testdf = testdf.drop(['sum'], axis=1) 

testdf.loc[:,['CladeA','CladeB','CladeC','CladeD','CladeE','CladeF','CladeG','CladeH', 'CladeI']].plot.line( color=colors, figsize=(5,5),linewidth=3,legend=None)


for key, value in colorsdict.items():
	print(key)
	print(value)

	#print(cidf[key])
	ci = stats.binom_conf_interval(cidf[key], trials, confidence_level=0.95)
	ci = pd.DataFrame(ci)
	plt.fill_between(cidf.index,ci.loc[0,:],ci.loc[1,:], color=value, alpha=.1)


plt.ylim([0,.7])
plt.gcf().subplots_adjust(bottom=0.15)
plt.margins(y=0,x=0)
plt.locator_params(axis="x", nbins=4)
plt.savefig('CladesbydateAmericaspercenttest.pdf', dpi=300,format='pdf',bbox_inches="tight")




