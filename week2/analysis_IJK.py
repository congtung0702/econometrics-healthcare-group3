"""
Assignment 2 - analyse voor opdracht I en J (en ondersteuning K).

Leest de output van crc_screening_run_assignment2.py (events.csv en
durations.csv) en berekent:
  I) gewonnen levensjaren door screening
  J) voorkomen CRC-gevallen en CRC-sterfgevallen
  
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DATA_DIR = Path(__file__).parent
EVENTS_FILE = DATA_DIR / "events.csv"
DURATIONS_FILE = DATA_DIR / "durations.csv"
OUTPUT_DIR = DATA_DIR / "output_IJK"
OUTPUT_DIR.mkdir(exist_ok=True)

events = pd.read_csv(EVENTS_FILE, index_col=0)
durations = pd.read_csv(DURATIONS_FILE, index_col=0)

# Totalen per tag per universe
ev = events.pivot_table(index="tag", columns="universe", values="number",
                        aggfunc="sum", fill_value=0)
du = durations.pivot_table(index="tag", columns="universe", values="number",
                           aggfunc="sum", fill_value=0)

NS, S = "no_screening", "screening"

print("=" * 60)
print("Totaal aantal events per universe")
print("=" * 60)
print(ev.to_string())

# ---------------------------------------------------------------------------
# I) Gewonnen levensjaren
# ---------------------------------------------------------------------------
ly_ns = du.loc["life", NS]
ly_s = du.loc["life", S]
ly_gained = ly_s - ly_ns
n = 100_000

print("\n" + "=" * 60)
print("I) Levensjaren")
print("=" * 60)
print(f"Levensjaren zonder screening : {ly_ns:,.0f}")
print(f"Levensjaren met screening    : {ly_s:,.0f}")
print(f"Gewonnen levensjaren         : {ly_gained:,.0f}")
print(f"Per persoon (n={n:,})       : {ly_gained / n:.3f} jaar "
      f"(= {ly_gained / n * 365.25:.0f} dagen)")

# ---------------------------------------------------------------------------
# J) Voorkomen gevallen en sterfgevallen
# ---------------------------------------------------------------------------
stages = ["1", "2", "3", "4"]
cl_tags = [f"crc_cl{s}" for s in stages]
sd_tags = [f"crc_sd{s}" for s in stages]


def get(tag, universe):
    return ev.loc[tag, universe] if tag in ev.index else 0


cases_ns = sum(get(t, NS) for t in cl_tags + sd_tags)
cl_s = sum(get(t, S) for t in cl_tags)
sd_s = sum(get(t, S) for t in sd_tags)
cases_s = cl_s + sd_s
deaths_ns = get("crc_death", NS)
deaths_s = get("crc_death", S)

cases_prev = cases_ns - cases_s
deaths_prev = deaths_ns - deaths_s

print("\n" + "=" * 60)
print("J) Voorkomen CRC-gevallen en -sterfte")
print("=" * 60)
print(f"CRC-gevallen zonder screening : {cases_ns:,}")
print(f"CRC-gevallen met screening    : {cases_s:,} "
      f"(klinisch {cl_s:,} + screen-detected {sd_s:,})")
print(f"Voorkomen gevallen            : {cases_prev:,} "
      f"({cases_prev / cases_ns:.1%})")
print()
print(f"CRC-sterfte zonder screening  : {deaths_ns:,}")
print(f"CRC-sterfte met screening     : {deaths_s:,}")
print(f"Voorkomen sterfgevallen       : {deaths_prev:,} "
      f"({deaths_prev / deaths_ns:.1%})")
print()
print(f"Onsets (laesies) zonder / met : {get('crc_onset', NS):,} / "
      f"{get('crc_onset', S):,}")
print(f"Adenomen verwijderd           : {get('crc_screening_adenoma', S):,}")
print(f"Kankers via screening gevonden: {get('crc_screening_cancer', S):,}")
print(f"OC-sterfte zonder / met       : {get('oc_death', NS):,} / "
      f"{get('oc_death', S):,}")
print(f"Gewonnen levensjaren per voorkomen CRC-sterfgeval: "
      f"{ly_gained / deaths_prev:.1f}")

# Stadiumverdeling
stage_table = pd.DataFrame({
    "no_screening": [get(f"crc_cl{s}", NS) + get(f"crc_sd{s}", NS) for s in stages],
    "screening": [get(f"crc_cl{s}", S) + get(f"crc_sd{s}", S) for s in stages],
}, index=[f"stadium {s}" for s in stages])
stage_pct = stage_table / stage_table.sum() * 100
print("\nStadiumverdeling (% van alle CRC-gevallen)")
print(stage_pct.round(1).to_string())

# Aantal screeningstests per ronde (nuttig voor K: belasting/kosten)
scr = events[(events.universe == S) & events.tag.str.startswith("crc_screening_")]
tests = scr.pivot_table(index="age", columns="tag", values="number",
                        aggfunc="sum", fill_value=0)
tests["totaal"] = tests.sum(axis=1)
print("\nScreeningsuitslagen per screeningsleeftijd")
print(tests.to_string())
print(f"Totaal aantal tests: {tests['totaal'].sum():,}")
print(f"Tests per gewonnen levensjaar: {tests['totaal'].sum() / ly_gained:.1f}")

# Samenvatting opslaan
summary = pd.DataFrame({
    "maat": ["levensjaren", "CRC-gevallen", "CRC-sterfte"],
    "no_screening": [ly_ns, cases_ns, deaths_ns],
    "screening": [ly_s, cases_s, deaths_s],
})
summary["verschil (screening - no_screening)"] = summary.screening - summary.no_screening
summary.to_csv(OUTPUT_DIR / "summary_IJ.csv", index=False)

# ---------------------------------------------------------------------------
# Figuren
# ---------------------------------------------------------------------------
# 1) CRC-gevallen en sterfte per leeftijd, beide universes
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
crc_case_tags = cl_tags + sd_tags
for universe, style in [(NS, "-"), (S, "--")]:
    sub = events[events.universe == universe]
    cases_age = sub[sub.tag.isin(crc_case_tags)].groupby("age").number.sum()
    deaths_age = sub[sub.tag == "crc_death"].groupby("age").number.sum()
    axes[0].plot(cases_age.index, cases_age.values, style, label=universe)
    axes[1].plot(deaths_age.index, deaths_age.values, style, label=universe)
axes[0].set_title("CRC-detections per ages")
axes[1].set_title("CRC-deaths per age")
for ax in axes:
    ax.set_xlabel("Age")
    ax.set_ylabel("Number")
    ax.legend()
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "J.png", dpi=150)

# 2) Totalen naast elkaar
fig, ax = plt.subplots(figsize=(6, 4.5))
x = range(2)
w = 0.35
ax.bar([i - w / 2 for i in x], [cases_ns, deaths_ns], w, label="no_screening")
ax.bar([i + w / 2 for i in x], [cases_s, deaths_s], w, label="screening")
ax.set_xticks(list(x))
ax.set_xticklabels(["CRC-gevallen", "CRC-sterfte"])
ax.set_ylabel("Aantal (per 100.000)")
ax.set_title("Effect van screening (50, 55, 60 jaar)")
ax.legend()
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "J_totals.png", dpi=150)

print(f"\nFiguren en samenvatting opgeslagen in: {OUTPUT_DIR}")
plt.show()
