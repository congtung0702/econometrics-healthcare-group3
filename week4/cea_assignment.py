# D,E : Imported packages
import pandas as pd
from itertools import product
import matplotlib.pyplot as plt

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

# D: Add all eight screening strategies
universes = [universe_screening, universe_no_screening]

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


# B: Duration tags for stage III and IV cancers
stage3_tags = tuple(
    f"crc_{diagnosis}3_{phase}"
    for diagnosis in ["cl", "sd"]
    for phase in ["init", "cont", "term", "term_oc"]
)

stage4_tags = tuple(
    f"crc_{diagnosis}4_{phase}"
    for diagnosis in ["cl", "sd"]
    for phase in ["init", "cont", "term", "term_oc"]
)

# B: Disutilities for each age group
dis3_young = 0.89 * (1 - 0.72 / 0.84)
dis4_young = 0.89 * (1 - 0.60 / 0.84)

# B: cea_data based on Quality of Life inputs.
cea_data = {
    "Screens": {
        "conditions": [],
        "events": {
            # A: Count all screens
            # C: Screening Costs and initial follow-up costs
            "crc_screening_cancer": [352.50 + 404., 0., 0., 1.],
            "crc_screening_negative": [352.50 + 0.084*421.95, 0., 0., 1.],
            "crc_screening_adenoma": [352.50 + 0.121*421.95, 0., 0., 1.]
        },
        "durations": {}
    },

    "Under_60": {
        "conditions": [("age", "<", 60)],
        "events": {},
        "durations": {
            "life": [0., 1., 0.89, 0.],
            stage3_tags: [0., 0., -dis3_young, 0.],
            stage4_tags: [0., 0., -dis4_young, 0.]
        }
    },

    "Age_60_plus": {
        "conditions": [("age", ">=", 60)],
        "events": {},
        "durations": {
            "life": [0., 1., 0.84, 0.],
            stage3_tags: [0., 0., -0.12, 0.],
            stage4_tags: [0., 0., -0.24, 0.]
        }
    },

    "Cancer_Costs" : {
        "conditions": [],
        "events": {},
        "durations": {
            # C: Initial care
            r"^crc_(cl|sd)[12]_init$":
                [3141.59 * 12, 0., 0., 0.],
            r"^crc_(cl|sd)[34]_init$":
                [8000.00 * 12, 0., 0., 0.],

            # C: Continuing care
            r"^crc_(cl|sd)[12]_cont$":
                [141.42 * 12, 0., 0., 0.],
            r"^crc_(cl|sd)[34]_cont$":
                [271.78 * 12, 0., 0., 0.],

            # C: Final year before death from other causes
            r"^crc_(cl|sd)[12]_term_oc$":
                [141.42 * 12, 0., 0., 0.],
            r"^crc_(cl|sd)[34]_term_oc$":
                [271.78 * 12, 0., 0., 0.],

            # C: Terminal care before CRC death
            r"^crc_(cl|sd)[1-4]_term$":
                [2109.80 * 12, 0., 0., 0.]
        }
    }
}

# C: Adjusted discount parameters
discount_rate = (0.03, 0.03, 0.03, 0.0)
discount_start = 2022

cea_result = cea.cost_effectiveness_analysis(
    cea_data=cea_data,
    cea_labels=cea_labels,
    durations=durations,
    events=events,
    discount_by="year",
    # C: set correct year
    discount_start=discount_start,
    discount_rate=discount_rate,
    reference="no_screening",
    # C: use regex notation
    use_regex=True
)

print(cea_result.cea_output)

# D: Obtain output off al 8 strategies
universe_table = pd.pivot_table(cea_result.cea_output,
aggfunc="sum", columns="label", values="number",
index="universe")

print(universe_table)

# E: Net Monetary Benefit
universe_table["NMB"] = (
    50000 * universe_table["QALY"] - universe_table["costs"]
)
print(universe_table[["NMB"]])

# E: Plot NMB for each strategy
ax = (universe_table["NMB"] / 1e6).sort_values().plot(
    kind="barh",
    figsize=(10, 5)
)

# E: Efficient frontier
efficient_strategies = cea.find_efficient(universe_table)
print(efficient_strategies)

# E: Cost-effectiveness scatter plot
cea.scatter_ce(
    universe_table,
    x_label="costs",
    y_label="QALY",
    tags="index",
    save_location="CEA_scatter.png"
)
cea.show()

ax.set_xlabel("Net Monetary Benefit (million)")
ax.set_ylabel("Screening Strategy")
ax.set_title("NMB at a Willingness-to-Pay of €50,000/QALY")

plt.tight_layout()
plt.savefig("NMB.png", dpi=300)
plt.show()
