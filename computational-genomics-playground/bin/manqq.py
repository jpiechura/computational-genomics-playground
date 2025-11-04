import sys, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

inp = Path(sys.argv[1])  # gwas_all.tsv
p_gws = float(sys.argv[2])  # 5e-8

df = pd.read_csv(inp, sep='\t')
df = df[['CHR','BP','SNP','A1','BETA','STAT','P']].dropna()
df = df[(df['P']>0) & (df['P']<=1)]

# Lambda (genomic inflation)
# Chi2 = qchisq(1 - P, df=1). For df=1, median of chi2 is ~0.4549364231
from scipy.stats import chi2
chisq = chi2.isf(df['P'], 1)
lambda_gc = np.median(chisq) / 0.4549364231
Path('lambda.txt').write_text(f"{lambda_gc:.4f}\n")

# BH-FDR q-values
p = df['P'].values
m = len(p)
order = np.argsort(p)
rank = np.empty_like(order)
rank[order] = np.arange(1, m+1)
q = p * m / rank
# monotone
q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
q_bh = np.empty_like(q_sorted)
q_bh[order] = q_sorted
df['Q_BH'] = q_bh

df.to_csv('gwas_all.with_q.tsv', sep='\t', index=False)

# Top hits by genome-wide threshold
tops = df[df['P'] <= p_gws].sort_values('P')
tops.to_csv('top_hits.tsv', sep='\t', index=False)

# Manhattan plot
df = df.sort_values(['CHR','BP'])
df['-LOG10P'] = -np.log10(df['P'])

# Build chromosome offsets for a continuous x-axis
chroms = df['CHR'].astype(int).unique()
offset = 0
ticks = []
ticklabels = []
xpos = []
dfs = []
for c in chroms:
    d = df[df['CHR'].astype(int)==c].copy()
    if d.empty: 
        continue
    d['x'] = d['BP'] + offset
    ticks.append(d['x'].median())
    ticklabels.append(str(c))
    offset = d['x'].max()
    dfs.append(d)
df2 = pd.concat(dfs)

plt.figure(figsize=(12,4))
plt.scatter(df2['x'], df2['-LOG10P'], s=2)
plt.axhline(-np.log10(p_gws), linestyle='--')
plt.ylabel('-log10(P)')
plt.xticks(ticks, ticklabels)
plt.xlabel('Chromosome')
plt.title('Manhattan plot')
plt.tight_layout()
plt.savefig('manhattan.png', dpi=200)
plt.close()

# QQ plot
obs = -np.log10(np.sort(df['P'].values))
exp = -np.log10((np.arange(1, len(obs)+1) / (len(obs)+1)))
mval = max(obs.max(), exp.max())
plt.figure(figsize=(4,4))
plt.scatter(exp, obs, s=3)
plt.plot([0, mval], [0, mval], linestyle='--')
plt.xlabel('Expected -log10(P)')
plt.ylabel('Observed -log10(P)')
plt.title('QQ plot')
plt.tight_layout()
plt.savefig('qq.png', dpi=200)
plt.close()