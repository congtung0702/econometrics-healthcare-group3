from pathlib import Path

import pandas as pd

events = pd.read_csv(Path(__file__).parent / "events.csv")

# A: number of positive and negative screening tests
screening_tests = events[events["tag"].str.startswith("crc_screening")]
print(screening_tests.groupby("tag")["number"].sum())
