import numpy as np
import pandas as pd

nphen = 5
shared = 0.25
specific = 0.15

R = 25
N = 5*10**4
M = 10**4
C = 'all_and_ind_overlaps_uniform_gg_random_ge_0.1'

shared_m = int(shared*M)
specific_m = int(specific*M)

gg = np.loadtxt('configs/%s_gg.txt'%C)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('../out/%s_%s_%s_RESULTS.txt'%(N,M,C),sep='\t')

corr_df = df.loc[df.Metric == 'Correlation']
corr_df.columns = ['|Correlation|']+list(df.columns[1:])
corr_df['|Correlation|'] = np.abs(corr_df['|Correlation|'])

powr_df = df.loc[df.Metric == 'Power']
powr_df.columns = ['Power']+list(df.columns[1:])

fig,axs = plt.subplots(ncols=2,figsize=(18,12),sharey=False)
for ax,d,m in zip(axs,[corr_df,powr_df],['|Correlation|','Power']):
    sns.boxplot(data=d,x='Heritability',y=m,hue='Source',ax=ax)
    ax.set_title('%s'%m)

plt.savefig('../out/%s_%s_%s_RESULTS.png'%(N,M,C))
