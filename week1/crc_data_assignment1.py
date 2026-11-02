"""
This files contains placeholder variables for the CRC process.
"""

# The individual onset age is drawn from a gamma distribution.
# This parameter holds the parameters for this distribution. A parameter
# set (i.e. 'shape' and 'scale') should be given.
# 1B:
age = {
    "shape": 30.0,
    "scale": 2.2
}

# Dictionary holding all the parameters determining the dwell times in each state.
# 1C:
dwell = {
    "adenoma": 10,
    "stage1" : 2.5,
    "stage2" : 2.5,
    "stage3" : 3.7,
    "stage4" : 1.5
}

# When a lesion transitions to a clinical state, it is cured
# (i.e. the individual will live with the disease until death from other
# causes) with a certain probability depending on the current age. This
# parameter holds PiecewiseLinear functions mapping ages to probabilities
# for each clinical stage.
cure = {
    "1": [
        (0, 0.905),
        (50, 0.881),
        (60, 0.860),
        (70, 0.753),
        (80, 0.682),
    ],
    "2": [
        (0, 0.705),
        (50, 0.658),
        (60, 0.601),
        (70, 0.493),
        (80, 0.424),
    ],
    "3": [
        (0, 0.561),
        (50, 0.556),
        (60, 0.539),
        (70, 0.452),
        (80, 0.376),
    ],
    "4": [
        (0, 0.088),
        (50, 0.064),
        (60, 0.050),
        (70, 0.032),
        (80, 0.023),
    ]
}

# If a cancerous lesion is not cured, an individual
# will die from the disease. This parameter holds the distributions from
# which the times to death are drawn. This parameter holds PiecewiseLinear
# functions mapping ages to probabilities for each clinical stage.
time_to_death = {
    "1": [
        (0.000, 0),
        (0.097, 1),
        (0.162, 2),
        (0.234, 3),
        (0.298, 4),
        (0.356, 5),
        (0.407, 6),
        (0.452, 7),
        (0.494, 8),
        (0.533, 9),
        (0.568, 10),
        (0.600, 11),
        (0.630, 12),
        (0.657, 13),
        (0.682, 14),
        (0.704, 15),
        (0.724, 16),
        (0.743, 17),
        (0.760, 18),
        (0.776, 19),
        (0.791, 20),
        (0.805, 21),
        (0.827, 22),
        (0.862, 23),
        (0.916, 24),
        (1.000, 25),
    ],
    "2": [
        (0.000, 0),
        (0.136, 1),
        (0.211, 2),
        (0.288, 3),
        (0.349, 4),
        (0.398, 5),
        (0.438, 6),
        (0.470, 7),
        (0.499, 8),
        (0.526, 9),
        (0.552, 10),
        (0.576, 11),
        (0.599, 12),
        (0.621, 13),
        (0.641, 14),
        (0.663, 15),
        (0.686, 16),
        (0.711, 17),
        (0.737, 18),
        (0.765, 19),
        (0.794, 20),
        (0.825, 21),
        (0.860, 22),
        (0.901, 23),
        (0.947, 24),
        (1.000, 25),
    ],
    "3": [
        (0.000, 0),
        (0.212, 1),
        (0.391, 2),
        (0.530, 3),
        (0.617, 4),
        (0.674, 5),
        (0.711, 6),
        (0.735, 7),
        (0.756, 8),
        (0.774, 9),
        (0.790, 10),
        (0.804, 11),
        (0.815, 12),
        (0.826, 13),
        (0.835, 14),
        (0.844, 15),
        (0.853, 16),
        (0.862, 17),
        (0.872, 18),
        (0.881, 19),
        (0.891, 20),
        (0.901, 21),
        (0.915, 22),
        (0.935, 23),
        (0.962, 24),
        (1.000, 25),
    ],
    "4": [
        (0.000, 0),
        (0.615, 1),
        (0.806, 2),
        (0.888, 3),
        (0.925, 4),
        (0.945, 5),
        (0.956, 6),
        (0.962, 7),
        (0.967, 8),
        (0.971, 9),
        (0.973, 10),
        (0.975, 11),
        (0.977, 12),
        (0.978, 13),
        (0.979, 14),
        (0.980, 15),
        (0.981, 16),
        (0.983, 17),
        (0.985, 18),
        (0.987, 19),
        (0.990, 20),
        (0.993, 21),
        (0.995, 22),
        (0.997, 23),
        (0.999, 24),
        (1.000, 25),
    ],
}
