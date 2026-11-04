# F: share of cancers detected by screening vs clinically, and stage distribution
from pathlib import Path

import pandas as pd

events = pd.read_csv(Path(__file__).parent / "eventsFG.csv", index_col=0)
scr = events[(events.universe == "screening") & events.tag.str.match(r"crc_(sd|cl)\d$")].copy()
scr["route"] = scr.tag.str[4:6].map({"sd": "screen-detected", "cl": "clinical"})
scr["stage"] = scr.tag.str[-1].astype(int)

counts = scr.pivot_table(index="stage", columns="route", values="number", aggfunc="sum", fill_value=0)
totals = counts.sum()

print("F: cancers per detection route (screening universe)")
table = counts.copy()
table["total"] = table.sum(axis=1)   # new column: both routes added, per stage
table.loc["all"] = table.sum()       # new row: all stages added, per column
print(table.to_string(), "\n")

print("Share of all detected cancers (%):")
print((100 * totals / totals.sum()).round(1).to_string(), "\n")

print("Stage distribution within each route (%):")
dist = 100 * counts / totals
dist["all (sd + cl)"] = 100 * counts.sum(axis=1) / totals.sum()
dist.loc["total"] = dist.sum()   
print(dist.round(1).to_string())
