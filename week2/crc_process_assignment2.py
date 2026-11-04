from itertools import product
from typing import Dict

import numpy as np

from miscore import Process, Individual
from miscore.utils import PiecewiseLinear

# When cancer is diagnosed by symptoms, this dictionary is used to translate
# each cancer state to a cancer stage
CANCER_STAGE = {                                                                                        # new
    "pcl4": "4",
    "pcl3b": "3",
    "pcl3a": "3",
    "pcl2b": "2",
    "pcl2a": "2",
    "pcl1b": "1",
    "pcl1a": "1",
}


class CRC(Process):
    """
    This is a reduced-form version of the panmodel CRC process that simulate
    colorectal lesions and cancer.

    :param name: The name of the process. Defaults to 'crc'.
    :param age: This parameter holds the parameters for the onset distribution.
        A parameter set (i.e. 'shape' and 'scale') should be given.
    :param dwell: Dictionary holding all the parameters determining the dwell
        times in each state.
    :param cure: When a lesion transitions to a clinical state, it is cured
        (i.e. the individual will live with the disease until death from other
        causes) with a certain probability depending on the current age. This
        parameter holds PiecewiseConstant functions mapping ages to probabilities.
    :param time_to_death: If a cancerous lesion is not cured, an individual
        will die from the disease. This parameter holds the distributions from
        which the times to death are drawn.
    """
    name = "crc"

    # log the important event tags
    event_tags = [
        "crc_onset",
        "crc_cl1",
        "crc_cl2",
        "crc_cl3",
        "crc_cl4",
        "crc_death",
    ]

    # log the important duration tags
    duration_tags = [
        f"crc_{sd_cl}{stage}_{phase}"
        for sd_cl, stage, phase
        in product(
            ["cl"],
            ["1", "2", "3", "4"],
            ["init", "cont", "term", "term_oc"]
        )
    ]

    def __init__(self,
                 age: Dict[str, float],
                 dwell: Dict[str, float],
                 cure: Dict[int, PiecewiseLinear.Data],
                 time_to_death: Dict[int, PiecewiseLinear.Data],
                 name: str = "crc"):
        self.name = name

        self.age = age
        self.dwell = dwell
        self.cure = {k: PiecewiseLinear(v) for k, v in cure.items()}
        self.time_to_death = {k: PiecewiseLinear(v) for k, v in time_to_death.items()}

        self.callbacks = {
            "__start__": self.schedule_onset,
            "crc_onset": self.onset,
            "crc_progress": self.progress,
            "crc_clinical": self.clinical,
            "__end__": self.log_durations
        }

    def properties(self, rng: np.random.Generator, n: int,
                   properties: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Method to define temporary properties for each individual."""
        # Calculate age of onset
        onset_age = rng.gamma(shape=self.age["shape"], scale=self.age["scale"], size=n)

        # Draw random numbers to use in the natural history of the lesions
        r_transition_1 = rng.random(n)
        r_transition_2 = rng.random(n)
        r_transition_3 = rng.random(n)
        r_dwell_precursor = rng.standard_exponential(n)
        r_cure = rng.random(n)
        r_time_to_death = rng.random(n)

        # B: Draw random numbers to use in sensitivity of the screening
        r_sensitivity_adenoma = rng.random(n)
        r_sensitivity_all_stage = rng.random(n)

        # B: Return the random numbers generated
        return {
            "crc_onset": onset_age,
            "crc_r_transition_1": r_transition_1,
            "crc_r_transition_2": r_transition_2,
            "crc_r_transition_3": r_transition_3,
            "crc_r_dwell": r_dwell_precursor,
            "crc_r_cure": r_cure,
            "crc_r_time_to_death": r_time_to_death,
            "crc_r_sensitivity_adenoma": r_sensitivity_adenoma,
            "crc_r_sensitivity_all_stage": r_sensitivity_all_stage,
        }

    @staticmethod
    def schedule_onset(individual: Individual):
        """Method to schedule disease at onset age. """

        # Initialize list that will hold the lesion
        individual.memory["crc_lesion"] = []                                                            # new

        # D: boolean that tells whether the lesion was removed (e.g. by screening)
        individual.memory["crc_lesion_removed"] = False

        # Schedule onset
        individual.add_to_queue_age(individual.properties("crc_onset"), "crc_onset")

    def onset(self, individual: Individual):
        """Method to create disease trajectory"""

        individual.log_event("crc_onset")

        # Retrieve natural history
        natural_history = self._natural_history(individual)
        state, dwell = next(natural_history)

        # Add event
        individual.add_to_queue(dwell, "crc_progress")

        # Save state machine
        individual.memory["crc_natural_history"] = natural_history

        # Store lesion in memory
        individual.memory["crc_lesion"] = state                                                         # new

    @staticmethod
    def progress(individual: Individual):
        """Method to simulate progression of the disease."""

        # D: a removed lesion no longer progresses
        if individual.memory["crc_lesion_removed"]:
            return

        # E: A clinically diagnosed cancer no longer progresses
        if "crc_clinical_information" in individual.memory:
            return

        # Retrieve natural history
        state, dwell = next(individual.memory["crc_natural_history"])

        # Schedule next progress event after the dwell time
        if dwell is not None:
            individual.add_to_queue(dwell, "crc_progress")
            # Store lesion in memory
            individual.memory["crc_lesion"] = state                                                     # new

        else:
            individual.add_to_queue(0., "crc_clinical", stage=state[-1])

    def clinical(self, individual: Individual, stage: str, screen_detected: bool = False):              # new
        """Method to determine clinical status of the lesions."""

        individual.memory["crc_clinical_information"] = {
            "age": individual.age,
            "stage": stage,
            "screen_detected": screen_detected                                                          # new
        }

        individual.log_event(f"crc_cl{stage}")

        # stop screening when clinically detected
        individual.add_to_queue(0., "crc_stop_screening")                                                  # new

        # Retrieve random numbers
        r_cure = individual.properties("crc_r_cure")
        cure_prob = self.cure[stage](individual.age)

        if r_cure < cure_prob:
            return

        # Retrieve random numbers
        r_time_to_death = individual.properties("crc_r_time_to_death")
        time_to_death = self.time_to_death[stage](r_time_to_death)

        individual.add_to_queue(time_to_death, "crc_death", terminate=True)

    def _natural_history(self, individual: Individual):
        """Method to set-up the natural history."""

        yield "adenoma", individual.properties("crc_r_dwell") * self.dwell["adenoma"]

        if individual.properties("crc_r_transition_1") < 0.14:
            yield "pcl1a", individual.properties("crc_r_dwell") * self.dwell["pcl1a"]
            yield "cl1", None
            return

        yield "pcl1b", individual.properties("crc_r_dwell") * self.dwell["pcl1b"]

        if individual.properties("crc_r_transition_2") < 0.33:
            yield "pcl2a", individual.properties("crc_r_dwell") * self.dwell["pcl2a"]
            yield "cl2", None
            return

        yield "pcl2b", individual.properties("crc_r_dwell") * self.dwell["pcl2b"]

        if individual.properties("crc_r_transition_3") < 0.49:
            yield "pcl3a", individual.properties("crc_r_dwell") * self.dwell["pcl3a"]
            yield "cl3", None
            return

        yield "pcl3b", individual.properties("crc_r_dwell") * self.dwell["pcl3b"]

        yield "pcl4", individual.properties("crc_r_dwell") * self.dwell["pcl4"]

        yield "cl4", None
        return

    @staticmethod
    def log_durations(individual: Individual, message: str):
        """Method to log the durations."""

        # Only log durations if CRC was clinical
        if "crc_clinical_information" not in individual.memory:
            return

        # Obtain the information on the clinical state of CRC of this person
        clinical = individual.memory["crc_clinical_information"]

        # Determine clinical stage, depending on whether or not CRC was
        # screen-detected
        if clinical["screen_detected"]:                                                                 # new
            stage = f"sd{clinical['stage']}"                                                            # new
        else:
            stage = f"cl{clinical['stage']}"

        # Obtain age at which CRC became clinical and age of death
        t0 = clinical["age"]
        t3 = individual.age

        # Calculate ages in between to be able to split the duration into init,
        # cont and term durations
        t2 = max(t0, t3 - 1.)
        t1 = min(t2, t0 + 1.)

        # Log durations
        individual.log_duration(f"crc_{stage}_init", t0, t1)
        individual.log_duration(f"crc_{stage}_cont", t1, t2)
        if message == "crc_death":
            individual.log_duration(f"crc_{stage}_term", t2, t3)
        else:
            individual.log_duration(f"crc_{stage}_term_oc", t2, t3)

