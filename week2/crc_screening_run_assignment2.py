from miscore import Model, processes, Universe
import crc_data_assignment2
import crc_screening_data_assignment2
from crc_process_assignment2 import CRC
from crc_screening_process_assignment2 import CRC_screening

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
    age=crc_data_assignment2.age,
    dwell=crc_data_assignment2.dwell,
    cure=crc_data_assignment2.cure,
    time_to_death=crc_data_assignment2.time_to_death
)

# CRC screening process
crc_screening = CRC_screening(
    name="crc_screening",
    sensitivity=crc_screening_data_assignment2.sensitivity,
    screening_ages=crc_screening_data_assignment2.screening_ages
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
screening = Universe(
    name="screening",
    processes=[
        birth,
        oc,
        crc,
        crc_screening
    ]
)

# Model
model = Model(
    universes=[
        no_screening, screening
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
