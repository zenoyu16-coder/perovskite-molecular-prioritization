"""Plot the completed CREST conformer-energy landscape of the matched isomers."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/quantum_pilot"
df=pd.read_csv(OUT/"crest_conformers.csv")
fig,ax=plt.subplots(figsize=(7.5,4.6))
for group,marker in (("sulfonyl_meta","o"),("sulfonyl_ortho","s")):
    subset=df[df.group.eq(group)]
    ax.scatter(subset.min_NO_A,subset.relative_kcal_mol,label=group.replace("_"," "),
               marker=marker,s=56,edgecolors="white",linewidths=.6,zorder=3)
ax.axvline(3.3,color="gray",alpha=.5,linestyle=":",linewidth=1.1)
ax.set(xlabel="Shortest intramolecular N···O distance (Å)",
       ylabel="Relative GFN2-xTB energy (kcal mol⁻¹)",
       title="Sulfonylbenzylammonium positional isomers")
ax.legend(frameon=False)
ax.grid(alpha=.18)
fig.tight_layout()
fig.savefig(OUT/"crest_isomer_conformations.png",dpi=200)
plt.close(fig)
print("Saved",OUT/"crest_isomer_conformations.png")
