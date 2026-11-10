import csv
import datetime
import numpy as np
from scipy.optimize import minimize

from miscore import Model, processes, Universe
from assignment3 import crc_data
from assignment3.crc_process import CRC
from assignment3 import crc_screening_data
from assignment3.crc_screening_process import CRC_screening
from miscore.tools.calibration import Target, Fit, MultinomialDeviance, PoissonDeviance, events_from_result, durations_from_result


def model_from_x(x):
    # Birth process
    birth = processes.Birth(
        name="birth",
        year=0
    )

    # OC process
    oc = processes.OC(
        name="oc",
        life_table=processes.oc.data.nl_2019.life_table
    )

    age_from_x = {
        "shape": x[0],
        "scale": x[1]
    }

    clinical_prob_from_x = {

    }

    # CRC process
    crc = CRC(
        name="crc",
        age=crc_data.age,
        dwell=crc_data.dwell,
        clinical_probability=crc_data.clinical_probability,
        cure=crc_data.cure,
        time_to_death=crc_data.time_to_death,
    )

    # CRC screening process
    crc_screening = CRC_screening(
        name="crc_screening",
        sensitivity=crc_screening_data.sensitivity,
        screening_ages=[]
    )

    screening = Universe(
        name="screening",
        processes=[birth, oc, crc, crc_screening]
    )

    model = Model(
        universes=[screening]
    )

    return model


target_incidence = Target(
    name="Incidence",
    x=[52, 57, 62, 67, 72, 77, 82, 93],
    obs=[922, 1908, 3345, 4789, 5774, 6149, 5577, 5966],
    n=[483074, 473140, 456751, 430442, 390066, 334394, 259603, 260922],
    result_to_sim=events_from_result(
        tags=[["crc_cl1"]],
        ages=[[], [], ...]
    ),
    result_to_m=durations_from_result(
        tags=[["life"]],
        ages=[[], [], ...]
    ),
    goodness_of_fit=PoissonDeviance(),
    weight=4.
)

target_stage = Target(
    name="Stage distribution",
    category_names=["Stage 1", "Stage 2", "Stage 3", "Stage 4"],
    x=[52, 57, 62, 67, 72, 77, 82, 93],
    obs=np.array([
        [149, 258, 256, 259],
        [306, 610, 504, 488],
        [554, 1060, 866, 865],
        [769, 1455, 1254, 1311],
        [957, 1816, 1566, 1435],
        [970, 1858, 1701, 1620],
        [817, 1707, 1545, 1508],
        [866, 1735, 1671, 1694],
    ]),
    n=[922, 1908, 3345, 4789, 5774, 6149, 5577, 5966],
    result_to_sim=events_from_result(
        tags=[["crc_cl1"]],
        ages=[[], [], ...]
    ),
    result_to_m=events_from_result(
        tags=[["crc_cl1"]],
        ages=[[], [], ...]
    ),
    goodness_of_fit=MultinomialDeviance()
)

fit = Fit(
    model_from_x=model_from_x,
    targets=[target_incidence, target_stage],
    n=1e5,
    seed=1234,
    event_ages=[50, 55, 60, 65, 70, 75, 80, 85],
    duration_ages=[50, 55, 60, 65, 70, 75, 80, 85],
)

# Your initial guess for the parameter set
x0 = [21.0, 3.9, 0.2, 0.22, 0.6]

# Only evaluate x0 -> Comment this out for question 2G
fit(x0)

Nfeval = 1
def callback_func(x):
    global Nfeval
    print(Nfeval)
    Nfeval += 1

# Calibrate the model -> Uncomment this for question 2G
# x0_2 = [21.0, 3.9, 0.2, 0.22, 0.6]
# print(datetime.datetime.now())
# minimize(
#     fun=fit,
#     x0=np.array(x0_2),
#     callback=callback_func,
#     method="Nelder-Mead",
#     options={
#         "adaptive": True,
#         "maxiter": 20,
#         "disp": True
#     }
# )
# print(datetime.datetime.now())


# Generate output - export the best found parameter set
with open(f"best.csv", "w+", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(
        ["fit", "x0", "x1", "x2", "x3", "x4"]
    )
    for i in fit.best:
        writer.writerow([i[0]] + [par for par in i[1]])

# Obtain the best found parameter set
best_result = fit.best[0][2]

# Generate a plot of the incidence target
target_incidence.plot(
    result=best_result,
    confidence_interval=target_incidence.confidence_interval(
        alpha=.05,
        distribution="poisson",
        method="normal"
    )
).savefig("incidence.png")

# Generate a plot of the stage distribution target
target_stage.plot(
    result=best_result,
    confidence_interval=target_stage.confidence_interval(
         alpha=.05,
         distribution="poisson",
         method="normal"
    )
).savefig("stage.png")

print("done")
