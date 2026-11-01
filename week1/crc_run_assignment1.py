from miscore import Model, processes, Universe
from assignment1 import crc_data_assignment1
from assignment1.crc_process_assignment1 import CRC

# Birth process
birth = processes.Birth(
    name="birth",
    year=1975
)

# OC process
oc = processes.OC(
    name="oc",
    life_table=processes.oc.data.nl_2019.life_table
)

# CRC process
crc = CRC(
    name="crc",
    age=crc_data_assignment1.age,
    dwell=crc_data_assignment1.dwell,
    cure=crc_data_assignment1.cure,
    time_to_death=crc_data_assignment1.time_to_death
)

# Universes
no_screening = Universe(
    name="no_screening",
    processes=[
        birth,
        oc,
        crc
    ]
)

# Model
model = Model(
    universes=[
        no_screening
    ]
)

# Run model
result = model.run(
    n=1e5,
    seed=1234,
    event_ages=range(0, 101, 1),
    event_years=range(1975, 2076, 1),
    duration_ages=range(0, 101, 1),
    duration_years=range(1975, 2076, 1),
    verbose=True
)

# Export
result.events.to_csv("events.csv")
result.durations.to_csv("durations.csv")
