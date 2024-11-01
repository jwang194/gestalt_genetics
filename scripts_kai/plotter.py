import numpy as np
import pandas as pd

nphen = 10
shared = 0.2
specific = 0.08

R = 25
N = 50000
M = 10000
C = 'traits_10_shared_0.2_specific_0.08_uniform_rg_0.05_random_re_-0.1'

shared_m = int(shared*M)
specific_m = int(specific*M)

gg = np.loadtxt('configs/%s_gg.txt'%C)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('/u/home/k/kaia/GESTALT/data/sim/results/%s_%s_%s_RESULTS.txt'%(N,M,C),sep='\t')

corr_df = df.loc[df.Metric == 'Correlation']
corr_df.columns = ['|Correlation|']+list(df.columns[1:])
corr_df['|Correlation|'] = np.abs(corr_df['|Correlation|'])

powr_df = df.loc[df.Metric == 'Power']
powr_df.columns = ['Power']+list(df.columns[1:])

fig,axs = plt.subplots(ncols=2,figsize=(18,12),sharey=False)
for ax,d,m in zip(axs,[corr_df,powr_df],['|Correlation|','Power']):
    sns.boxplot(data=d,x='Heritability',y=m,hue='Source',ax=ax)
    ax.set_title('%s'%m)

plt.savefig('/u/home/k/kaia/GESTALT/data/sim/results/%s_%s_%s_RESULTS.png'%(N,M,C))
