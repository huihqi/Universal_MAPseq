import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
M='matrix'
cols = np.r_[1:8, 9:16]

sub = M[:, cols]                      # selected columns
row_keep = (sub > 0).any(axis=1)      # at least one > 0 in those columns
merged = sub[row_keep, :] 


#1) Calculate P(B|A)
X = np.asarray(merged)
Xb = (X > 0).astype(np.int64)

nA = Xb.sum(axis=0)
co = Xb.T @ Xb

P = co / nA[:, None].astype(float)
P[nA == 0, :] = np.nan


if hasattr(merged, "columns"):
    names = list(merged.columns)
else:
    names = [f"R{i+1}" for i in range(Xb.shape[1])]

P_df = pd.DataFrame(P, index=names, columns=names)
P_df.to_csv("P_B_given_A.csv", index=True)

#  2) Save heatmap
fig, ax = plt.subplots(figsize=(10, 8))
rgb = np.loadtxt('/Users/huihqi/Downloads/Vocalization/parula_256.csv', delimiter=",")  # shape (256, 3)
parula = ListedColormap(rgb, name="parula")
parula.set_bad(color="white")
im = ax.imshow(P_df.to_numpy(), vmin=0, vmax=1, aspect="auto", cmap=parula)

ax.set_title("P(B|A) heatmap")
ax.set_xlabel("B (columns)")
ax.set_ylabel("A (rows)")

ax.set_xticks(np.arange(len(names)))
ax.set_yticks(np.arange(len(names)))
ax.set_xticklabels(names, rotation=90)
ax.set_yticklabels(names)

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("P(B|A)")

fig.tight_layout()
fig.savefig("P_B_given_A_heatmap_Parrot_2L_wo_NLC.eps", dpi=300, bbox_inches="tight")
plt.show()