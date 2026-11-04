from pathlib import Path

import pandas as pd

events = pd.read_csv(Path(__file__).parent / "events.csv")

# A: number of positive and negative screening tests
screening_tests = events[events["tag"].str.startswith("crc_screening")]
print(screening_tests.groupby("tag")["number"].sum())

# D: clinical cancers and CRC deaths per universe (removed adenomas no longer progress)
crc_events = events[events["tag"].isin(["crc_cl1", "crc_cl2", "crc_cl3", "crc_cl4", "crc_death"])]
print(crc_events.pivot_table(index="tag", columns="universe", values="number", aggfunc="sum"))
