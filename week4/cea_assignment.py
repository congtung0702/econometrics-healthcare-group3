from itertools import product

import matplotlib.pyplot as plt
import pandas as pd

from miscore import Model, processes, Universe
import crc_data, crc_screening_data
from miscore.tools.cea import cea
from cea_processes import CRC, CRC_screening

birth = processes.Birth(
    name="birth",
    year=1972
)

# OC process
oc = processes.OC(
    name="oc",
    life_table=processes.oc.data.nl_2019.life_table
)

# CRC process
crc = CRC(
    name="crc",
    age=crc_data.age,
    dwell=crc_data.dwell,
    cure=crc_data.cure,
    time_to_death=crc_data.time_to_death,
    clinical_probability=crc_data.clinical_probability
)

crc_screening = CRC_screening(
    name="crc_screening",
    sensitivity=crc_screening_data.sensitivity,
    screening_ages=[50, 52, 54]
)

universe_screening = Universe(
    name="screening",
    processes=[birth, oc, crc, crc_screening]
)

universe_no_screening = Universe(
    name="no_screening",
    processes=[birth, oc, crc]
)

# Define which strategies are included in our CEA
universes = [
    universe_screening,
    universe_no_screening
]

# D: Generate all eight screening strategies
for i, (start, interval, stop) in enumerate(
    product([50, 55], [2, 3], [60, 65]), start=1
):
    screen = CRC_screening(
        name=f"screening_{i}",
        sensitivity=crc_screening_data.sensitivity,
        screening_ages=list(range(start, stop, interval))
    )

    universes.append(
        Universe(
            name=f"universe_screening_{i}",
            processes=[birth, oc, crc, screen]
        )
    )

# Pass the list of universes to our Model
model = Model(universes)

# Run our model, tallying events and durations from age 40 onwards
result = model.run(
    n=1e5,
    seed=42195,
    event_years=range(2022, 2072, 1),
    duration_years=range(2022, 2072, 1),
    duration_ages=range(40, 101, 1),
    verbose=True
)

result.events.to_csv("cea_events.csv")
result.durations.to_csv("cea_durations.csv")

events = result.events
durations = result.durations

cea_labels = [
    "costs",
    "LY",
    "QALY",
    "screens"
]

# B: Disutilities for individuals younger than 60 (same relative reduction)
dis3_young = 0.89 * (1 - 0.72 / 0.84)
dis4_young = 0.89 * (1 - 0.60 / 0.84)

# B: Regex patterns for all duration tags of stage III and IV
# B: (anchored with ^...$ so that e.g. "term" does not also match "term_oc")
STAGE3 = r"^crc_(cl|sd)3_(init|cont|term|term_oc)$"
STAGE4 = r"^crc_(cl|sd)4_(init|cont|term|term_oc)$"

cea_data = {
    "Screens": {
        "conditions": [],
        "events": {
            # A: Count all screens
            # C: Screening costs (invitation 2.50 + screen 350) and follow-up costs
            r"^crc_screening_cancer$": [352.50 + 404., 0., 0., 1.],
            r"^crc_screening_negative$": [352.50 + 0.084 * 421.95, 0., 0., 1.],
            r"^crc_screening_adenoma$": [352.50 + 0.121 * 421.95, 0., 0., 1.]
        },
        "durations": {}
    },

    # B: Norm utility and disutilities below age 60
    "Under_60": {
        "conditions": [("age", "<", 60)],
        "events": {},
        "durations": {
            r"^life$": [0., 1., 0.89, 0.],
            STAGE3: [0., 0., -dis3_young, 0.],
            STAGE4: [0., 0., -dis4_young, 0.]
        }
    },

    # B: Norm utility and disutilities from age 60
    "Age_60_plus": {
        "conditions": [("age", ">=", 60)],
        "events": {},
        "durations": {
            r"^life$": [0., 1., 0.84, 0.],
            STAGE3: [0., 0., -0.12, 0.],
            STAGE4: [0., 0., -0.24, 0.]
        }
    },

    # C: Cancer treatment costs (monthly costs * 12, durations are in years)
    "Cancer_Costs": {
        "conditions": [],
        "events": {},
        "durations": {
            # C: Initial care
            r"^crc_(cl|sd)[12]_init$": [3141.59 * 12, 0., 0., 0.],
            r"^crc_(cl|sd)[34]_init$": [8000.00 * 12, 0., 0., 0.],
            # C: Continuing care
            r"^crc_(cl|sd)[12]_cont$": [141.42 * 12, 0., 0., 0.],
            r"^crc_(cl|sd)[34]_cont$": [271.78 * 12, 0., 0., 0.],
            # C: Final year before death from other causes -> continuing care
            r"^crc_(cl|sd)[12]_term_oc$": [141.42 * 12, 0., 0., 0.],
            r"^crc_(cl|sd)[34]_term_oc$": [271.78 * 12, 0., 0., 0.],
            # C: Terminal care before CRC death
            r"^crc_(cl|sd)[1-4]_term$": [2109.80 * 12, 0., 0., 0.]
        }
    }
}

# C: 3% discounting for costs, LY and QALY (not for screens), from 2022
discount_rate = (0.03, 0.03, 0.03, 0.0)
discount_start = 2022

cea_result = cea.cost_effectiveness_analysis(
    cea_data=cea_data,
    cea_labels=cea_labels,
    durations=durations,
    events=events,
    discount_by="year",
    discount_start=discount_start,
    discount_rate=discount_rate,
    reference="no_screening",
    use_regex=True
)

print(cea_result.cea_output)

# D: Obtain output of all screening strategies
universe_table = pd.pivot_table(
    cea_result.cea_output,
    aggfunc="sum",
    columns="label",
    values="number",
    index="universe"
)

# C: ICERs of the original strategy versus no screening
print("ICER per QALY:", universe_table.loc["screening", "costs"] / universe_table.loc["screening", "QALY"])
print("ICER per LY:  ", universe_table.loc["screening", "costs"] / universe_table.loc["screening", "LY"])

# E: Net Monetary Benefit at a WTP of 50,000 per QALY
universe_table["NMB"] = 50000 * universe_table["QALY"] - universe_table["costs"]

print(universe_table)
universe_table.to_csv("universe_table.csv")

# E: Plot NMB for each strategy
ax = (universe_table["NMB"] / 1e6).sort_values().plot(
    kind="barh",
    figsize=(10, 5)
)
ax.set_xlabel("Net Monetary Benefit (€ million)")
ax.set_ylabel("Screening Strategy")
ax.set_title("NMB at a Willingness-to-Pay of €50,000/QALY")
plt.tight_layout()
plt.savefig("NMB.png", dpi=300, bbox_inches="tight")
plt.show()

# E: Cost-effectiveness frontier
tab = universe_table[["costs", "QALY"]].copy()
# The reference strategy (no screening) has 0 incremental costs and QALYs
tab.loc["no_screening"] = [0.0, 0.0]
tab = tab.sort_values("costs")

status = {}

# 1) Strong dominance: another strategy is cheaper (or equal) AND more effective
for name, row in tab.iterrows():
    if ((tab["costs"] <= row["costs"]) & (tab["QALY"] > row["QALY"])).any():
        status[name] = "dominated"

# 2) Extended dominance: ICER higher than that of the next strategy on the frontier
front = tab.drop(index=list(status)).copy()
while True:
    front = front.sort_values("costs")
    icer = front["costs"].diff() / front["QALY"].diff()
    bad = [front.index[i] for i in range(1, len(front) - 1)
           if icer.iloc[i] > icer.iloc[i + 1]]
    if not bad:
        break
    status[bad[0]] = "extendedly dominated"
    front = front.drop(index=bad[0])

front["ICER"] = front["costs"].diff() / front["QALY"].diff()
print(front)
print(status)

fig, ax = plt.subplots(figsize=(7, 5))
for name, row in tab.iterrows():
    colour = "C0" if name in front.index else "grey"
    ax.scatter(row["QALY"], row["costs"] / 1e6, color=colour, zorder=3)
    ax.annotate(name, (row["QALY"], row["costs"] / 1e6),
                textcoords="offset points", xytext=(5, -10), fontsize=8)
ax.plot(front["QALY"], front["costs"] / 1e6, "-", color="C0",
        label="Efficient frontier")
ax.set_xlabel("QALYs gained (discounted, per 100,000)")
ax.set_ylabel("Incremental costs (€ million, discounted)")
ax.set_title("Cost-effectiveness frontier")
ax.legend()
fig.tight_layout()
fig.savefig("frontier.png", dpi=300)
plt.show()