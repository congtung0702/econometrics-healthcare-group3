from itertools import product
from typing import Dict

from math import ceil
import numpy as np

from miscore import MAX_AGE
from miscore import Individual, Process
from miscore.utils import PiecewiseLinear


class CRC(Process):
    """
    This is a reduced-form version of the MISCAN CRC process that simulates
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

    event_tags = [
        "crc_pcl1a",
        "crc_pcl1b",
        "crc_pcl2a",
        "crc_pcl2b",
        "crc_pcl3a",
        "crc_pcl3b",
        "crc_pcl4",
        "crc_cl1",
        "crc_cl2",
        "crc_cl3",
        "crc_cl4",
        "crc_sd1",
        "crc_sd2",
        "crc_sd3",
        "crc_sd4",
        "crc_death",
    ]

    duration_tags = [
        f"crc_{sd_cl}{stage}_{phase}"
        for sd_cl, stage, phase
        in product(
            ["cl", "sd"],
            ["1", "2", "3", "4"],
            ["init", "cont", "term", "term_oc"]
        )
    ]

    def __init__(self,
                 age: Dict[str, float],
                 dwell: Dict[str, float],
                 cure: Dict[int, PiecewiseLinear.Data],
                 time_to_death: Dict[int, PiecewiseLinear.Data],
                 clinical_probability: Dict[int, int],
                 name: str = "crc"):
        self.name = name

        self.age = age
        self.dwell = dwell
        self.cure = {k: PiecewiseLinear(v) for k, v in cure.items()}
        self.time_to_death = {k: PiecewiseLinear(v) for k, v in time_to_death.items()}
        self.clinical_probability = clinical_probability

        self.callbacks = {
            "__start__": self.schedule_onset,
            "crc_onset": self.onset,
            "crc_progress": self.progress,
            "crc_clinical": self.clinical,
            "__end__": self.log_durations
        }

    def properties_tmp(self, rng: np.random.Generator, n: int,
                       properties: Dict[str, np.ndarray],
                       properties_tmp: Dict[str, np.ndarray]
                       ) -> Dict[str, np.ndarray]:
        """Method to define temporary properties for each individual."""
        # Calculate age of onset
        onset_age = rng.gamma(shape=self.age["shape"], scale=self.age["scale"], size=n)

        # 70% of individuals does not get an adenoma
        no_adenoma = rng.random(n) < 0.7
        onset_age[no_adenoma] = 100

        # Draw random numbers to use in the natural history of the lesions
        r_transition_1 = rng.random(n)
        r_transition_2 = rng.random(n)
        r_transition_3 = rng.random(n)
        r_cure = rng.random(n)
        r_time_to_death = rng.random(n)

        r_states = rng.standard_exponential((n, 8))

        return {
            "crc_onset": onset_age,
            "crc_r_transition_1": r_transition_1,
            "crc_r_transition_2": r_transition_2,
            "crc_r_transition_3": r_transition_3,
            "crc_r_cure": r_cure,
            "crc_r_time_to_death": r_time_to_death,
            "crc_r_states": r_states
        }

    def schedule_onset(self, simulation: Individual):
        """Method to schedule disease at onset age. """

        # Initialize list that will hold the lesion
        simulation.memory["crc_lesion"] = []
        simulation.memory["crc_lesion_removed"] = False

        # Retrieve onset age
        onset_age = simulation.properties_tmp("crc_onset")

        # Schedule onset
        if onset_age < MAX_AGE:
            simulation.add_to_queue_age(onset_age, "crc_onset")

    def onset(self, individual: Individual):
        """Method to create disease trajectory"""

        natural_history = self._natural_history(individual)
        state, dwell = next(natural_history)

        individual.add_to_queue(dwell, "crc_progress")
        individual.memory["crc_progress_age"] = individual.age + (dwell if dwell is not None else 0)

        # Store lesion in memory
        individual.memory["crc_lesion"] = {
            "state": state,
            "natural_history": natural_history,
        }

    @staticmethod
    def progress(individual: Individual):
        """Method to simulate progression of the disease."""
        lesion = individual.memory["crc_lesion"]

        # Do nothing if the lesion was removed
        if "removed" in lesion:
            return
        if "crc_clinical" in individual.memory:
            return

        state, dwell = next(lesion["natural_history"])

        lesion["state"] = state

        individual.memory["crc_progress_age"] = individual.age + (dwell if dwell is not None else 0)

        # Schedule next progress event after the dwell time
        if dwell is not None:
            individual.add_to_queue(dwell, "crc_progress")
        else:
            individual.add_to_queue(0., "crc_clinical", stage=state[-1])

    def clinical(self, individual: Individual, stage: str,
                 screen_detected: bool = False):
        """Method to determine clinical status of the lesions."""
        if "removed" in individual.memory["crc_lesion"]:
            return

        if screen_detected:
            individual.log_event(f"crc_sd{stage}")
        else:
            individual.log_event(f"crc_cl{stage}")

        individual.memory["crc_clinical"] = {
            "age": individual.age,
            "stage": stage,
            "screen_detected": screen_detected
        }

        individual.add_to_queue(0., "crc_stop_screening")

        # Retrieve random numbers
        r_cure = individual.properties_tmp("crc_r_cure")
        cure_prob = self.cure[stage](individual.age)

        if r_cure < cure_prob:
            return

        # Retrieve random numbers
        r_time_to_death = individual.properties_tmp("crc_r_time_to_death")
        time_to_death = self.time_to_death[stage](r_time_to_death)
        if screen_detected:
            time_to_death += individual.memory["crc_progress_age"] - individual.age
            if stage != '4':
                while True:
                    state, dwell = next(individual.memory["crc_lesion"]['natural_history'], (None, None))
                    if dwell is None:
                        break
                    time_to_death += dwell

        individual.add_to_queue(time_to_death, "crc_death", terminate=True)

    def _natural_history(self, individual: Individual):

        yield "adenoma", individual.properties_tmp("crc_r_states")[0] * self.dwell["adenoma"]

        if individual.properties_tmp("crc_r_transition_1") < self.clinical_probability[1]:
            yield "pcl1a", individual.properties_tmp("crc_r_states")[1] * self.dwell["pcl1a"]
            yield "cl1", None
            return

        yield "pcl1b", individual.properties_tmp("crc_r_states")[2] * self.dwell["pcl1b"]

        if individual.properties_tmp("crc_r_transition_2") < self.clinical_probability[2]:
            yield "pcl2a", individual.properties_tmp("crc_r_states")[3] * self.dwell["pcl2a"]
            yield "cl2", None
            return

        yield "pcl2b", individual.properties_tmp("crc_r_states")[4] * self.dwell["pcl2b"]

        if individual.properties_tmp("crc_r_transition_3") < self.clinical_probability[3]:
            yield "pcl3a", individual.properties_tmp("crc_r_states")[5] * self.dwell["pcl3a"]
            yield "cl3", None
            return

        yield "pcl3b", individual.properties_tmp("crc_r_states")[6] * self.dwell["pcl3b"]

        yield "pcl4", individual.properties_tmp("crc_r_states")[7] * self.dwell["pcl4"]

        yield "cl4", None
        return

    @staticmethod
    def log_durations(individual: Individual, message: str):
        """Method to log the durations."""
        # Only log durations if CRC was clinical
        if "crc_clinical" not in individual.memory:
            return

        # Obtain the information on the clinical state of CRC of this person
        clinical = individual.memory["crc_clinical"]

        # Determine clinical stage, depending on whether or not CRC was
        # screen-detected
        if clinical["screen_detected"]:
            stage = f"sd{clinical['stage']}"
        else:
            stage = f"cl{clinical['stage']}"

        # Obtain age at which CRC became clinical and age of death
        t0 = clinical["age"]
        t3 = individual.age

        # Death age should not be before clinical age
        assert t3 >= t0

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


import random
from typing import Any, Dict, List

from miscore import Process, Individual

# All preclinical colorectal cancer states
PRECLINICAL_STATES = {
    "pcl4",
    "pcl3b",
    "pcl3a",
    "pcl2b",
    "pcl2a",
    "pcl1b",
    "pcl1a",
}

# When cancer is diagnosed by symptoms, this dictionary is used to translate
# each cancer state to a cancer stage
CANCER_STAGE = {
    "pcl4": "4",
    "pcl3b": "3",
    "pcl3a": "3",
    "pcl2b": "2",
    "pcl2a": "2",
    "pcl1b": "1",
    "pcl1a": "1",
}


class CRC_screening(Process):
    name = "crc_screening"

    event_tags = [
        "crc_screening_first_invitation",
        "crc_screening_invitation",
    ]
    event_tags += [
        "crc_screening_positive",
        "crc_screening_negative",
        "crc_screening_adenoma",                                                                        # 4E
        "crc_screening_cancer",
    ]

    def __init__(self,
                 sensitivity: Dict[str, Any],
                 screening_ages: List[int],
                 name: str = "crc_screening"):
        self.name = name

        self.sensitivity = sensitivity
        self.screening_ages = screening_ages

        self.callbacks = {
            "__start__": self.schedule_screening,
            "crc_screening_invitation": self.crc_screen,
            "crc_stop_screening": self.stop_screening
        }

    def properties_tmp(self, rng: np.random.Generator, n: int,
                       properties: Dict[str, np.ndarray],
                       properties_tmp: Dict[str, np.ndarray]
                       ) -> Dict[str, np.ndarray]:
        return {
            "crc_r_sensitivity": rng.random((n, 20))
        }
    @staticmethod
    def stop_screening(individual: Individual):
        individual.memory["crc_screening_stop"] = True

    def schedule_screening(self, individual: Individual):
        """Method to schedule the screening invitations."""
        individual.memory["screen_no"] = 0
        for scr_age in self.screening_ages:
            individual.add_to_queue_age(scr_age, "crc_screening_invitation")

    def crc_screen(self, individual: Individual):
        """Method to process participation and outcomes of invited screening round."""
        # Do nothing if screening was stopped
        if "crc_screening_stop" in individual.memory:
            return

        if individual.memory["crc_lesion"]:
            lesion = individual.memory["crc_lesion"]
            if individual.properties_tmp("crc_r_sensitivity")[individual.memory["screen_no"]] < self.sensitivity[lesion["state"]]:
                if lesion["state"] == "adenoma":
                    result = "adenoma"
                    individual.memory["crc_lesion_removed"] = True
                else:
                    result = "cancer"
                    individual.add_to_queue(0., "crc_clinical",
                                         stage=CANCER_STAGE[individual.memory["crc_lesion"]["state"]],
                                         screen_detected=True)
            else:
                result = "negative"
        else:
            result = "negative"

        individual.memory["screen_no"] += 1
        individual.log_event(f"crc_screening_{result}")
