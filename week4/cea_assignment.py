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

cea_data = {
    "Example": {
        "conditions": [
        ],
        "events": {
            "crc_screening_cancer": [0., 0., 0., 1.]
        },
        "durations": {
            "life": [0., 1., 1., 0.]
        }
    }
}

discount_rate = (0.05, 0.0, 0.0, 0.0)

discount_start = 0

cea_result = cea.cost_effectiveness_analysis(
    cea_data=cea_data,
    cea_labels=cea_labels,
    durations=durations,
    events=events,
    discount_by="year",
    discount_start=0,
    discount_rate=discount_rate,
    reference="no_screening"
)

print(cea_result.cea_output)