# G: number of screen-detected and clinically detected cancers against age
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

here = Path(__file__).parent
events = pd.read_csv(here / "eventsFG.csv", index_col=0)
scr = events[(events.universe == "screening") & events.tag.str.match(r"crc_(sd|cl)\d$")]

by_age = scr.assign(route=scr.tag.str[4:6]).pivot_table(
    index="age", columns="route", values="number", aggfunc="sum", fill_value=0
).reindex(range(0, 101), fill_value=0)

print("G: peak ages")
print("clinical:", by_age.cl.idxmax(), "(", by_age.cl.max(), ")")
print("screen-detected:", by_age.sd.idxmax(), "(", by_age.sd.max(), ")")
print("screen-detected ages with cases:", list(by_age.index[by_age.sd > 0]))

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(by_age.index, by_age.cl, label="Clinically detected")
ax.plot(by_age.index, by_age.sd, label="Screen-detected")
ax.set(xlabel="Age", ylabel="Number of cancers", title="Detected cancers by age (screening universe)")
ax.legend()
fig.tight_layout()
fig.savefig(here / "G_detected_by_age.png", dpi=150)
