import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base = Path('/mnt/data/hw4_unpacked/hw4/hw4_2')
v = pd.read_csv(base/'values.csv')
e = pd.read_csv(base/'relerr.csv')

# sanitize NaNs/infs for plotting
for df in (v,e):
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
        df[c].replace([np.inf,-np.inf], np.nan, inplace=True)

# 1 values plots per precision
specs = [
    ('float', ['direct_f','expand_f','horner_f']),
    ('double', ['direct_d','expand_d','horner_d']),
    ('quad', ['direct_q','expand_q','horner_q']),
]
for name, cols in specs:
    plt.figure(figsize=(8,5.2))
    for c in cols:
        plt.plot(v['x'], v[c], label=c, linewidth=1.5)
    plt.xlim(0.7,1.3)
    plt.xlabel('x')
    plt.ylabel('value')
    plt.title(f'(x-1)^10 on [0.7, 1.3] ({name})')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(base/f'values_{name}.png', dpi=180)
    plt.close()

# 2 relative error plots (cap huge values for readability? use raw with ylim)
for name, cols in [('float',['rel_direct_f','rel_expand_f','rel_horner_f']),('double',['rel_direct_d','rel_expand_d','rel_horner_d']),('quad',['rel_direct_q','rel_expand_q','rel_horner_q'])]:
    plt.figure(figsize=(8,5.2))
    for c in cols:
        y = e[c].to_numpy(dtype=float)
        y = np.where(y<=0, np.nan, y)
        plt.semilogy(e['x'], y, label=c, linewidth=1.2)
    plt.xlim(0.7,1.3)
    plt.xlabel('x')
    plt.ylabel('relative error')
    plt.title(f'Relative error on [0.7, 1.3] ({name})')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(base/f'relerr_{name}.png', dpi=180)
    plt.close()

# 3 focused comparison using double/quad to show expansion failure near x=1
mask = (v['x']>=0.97) & (v['x']<=1.03)
plt.figure(figsize=(8,5.2))
for c in ['direct_d','expand_d','horner_d','direct_q']:
    plt.plot(v.loc[mask,'x'], v.loc[mask,c], label=c, linewidth=1.5)
plt.xlabel('x')
plt.ylabel('value')
plt.title('Zoom near x=1 (double and quad reference)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(base/'values_zoom_double_quad.png', dpi=180)
plt.close()
