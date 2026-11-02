from typing import Dict

import numpy as np

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

    # 1F, 2D: log the important event tags
    event_tags = [
        "crc_onset",
        "crc_cl1",
        "crc_cl2",
        "crc_cl3",
        "crc_cl4",
        "crc_death"
    ]

    # 1F: log the important duration tags
    duration_tags = [
        "crc_cl1",
        "crc_cl2",
        "crc_cl3",
        "crc_cl4",
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

        # 1B: Draw the age of onset for each individual
        onset_age = rng.gamma(self.age["shape"], self.age["scale"], size=n)

        # 1C: Draw transition probabilities
        trans_prob_stage1 = rng.uniform(0, 1, size=n)
        trans_prob_stage2 = rng.uniform(0, 1, size=n)
        trans_prob_stage3 = rng.uniform(0, 1, size=n)

        # 1D: Draw dwell time variation
        dwell_variation_adenoma = rng.exponential(1, size=n)
        dwell_variation_stage1 = rng.exponential(1, size=n)
        dwell_variation_stage2 = rng.exponential(1, size=n)
        dwell_variation_stage3 = rng.exponential(1, size=n)
        dwell_variation_stage4 = rng.exponential(1, size=n)

        # 2B: Draw cure probabilities
        cure_prob = rng.uniform(0, 1, size=n)

        # 2C: Draw time of death probabilities
        time_till_death = rng.uniform(0, 1, size=n)

        # 1B-1D, 2B-2C: Return the calculated properties
        return {
            "onset_age": onset_age,
            "trans_prob_stage1": trans_prob_stage1,
            "trans_prob_stage2": trans_prob_stage2,
            "trans_prob_stage3": trans_prob_stage3,
            "dwell_variation_adenoma": dwell_variation_adenoma,
            "dwell_variation_stage1": dwell_variation_stage1,
            "dwell_variation_stage2": dwell_variation_stage2,
            "dwell_variation_stage3": dwell_variation_stage3,
            "dwell_variation_stage4": dwell_variation_stage4,
            "cure_prob": cure_prob,
            "time_till_death": time_till_death,
        }

    @staticmethod
    def schedule_onset(individual: Individual):
        # 1B: Schedule onset
        individual.add_to_queue_age(individual.properties("onset_age"), "crc_onset")

    def onset(self, individual: Individual):
        """Method to create disease trajectory"""

        # Log event and save to memory
        individual.log_event("crc_onset")
        individual.memory["crc_onset"] = {
            "age": individual.age,
            "stage": "preclinical"
        }

        # Retrieve natural history
        natural_history = self._natural_history(individual)
        state, dwell = next(natural_history)

        # Add event
        individual.add_to_queue(dwell, "crc_progress")

        # Save state machine
        individual.memory["crc_natural_history"] = natural_history

    def progress(self, individual: Individual):
        """Method to simulate progression of the disease."""

        # 1E: Get natural history with state and dwell time
        natural_history = individual.memory["crc_natural_history"]
        state, dwell = next(natural_history)

        # 1E: Clinical next stage
        if state in ["crc_cl1", "crc_cl2", "crc_cl3", "crc_cl4"]:
            self.clinical(individual, state[-1])

        else:
            individual.add_to_queue(dwell, "crc_progress")

    def clinical(self, individual: Individual, stage: str):
        """Method to determine clinical status of the lesions."""

        # 1F: Log event and save to memory
        individual.log_event(f"crc_cl{stage}")
        individual.memory["crc_information"] = {
            "age": individual.age,
            "stage": stage
        }

        # 2B: Check whether one is cured or not
        cure_prob = self.cure[stage](individual.age)
        if individual.properties("cure_prob") < cure_prob:
            return

        # 2C: Determine the time of death
        time_till_death_years = self.time_to_death[stage](individual.properties("time_till_death"))

        # 2D: Schedule a terminal event
        individual.add_to_queue(
            delta=time_till_death_years,
            message="crc_death",
            terminate=True
        )

    def _natural_history(self, individual: Individual):
        """Method to set up the natural history."""

        # 1C-1D: Describes the natural history of CRC disease process.
        yield "crc_adenoma", self.dwell["adenoma"] * individual.properties("dwell_variation_adenoma")

        if individual.properties("trans_prob_stage1") < 0.14:
            yield "crc_pcl1a", self.dwell["stage1"] * individual.properties("dwell_variation_stage1")
            yield "crc_cl1", None
            return
        else:
            yield "crc_pcl1b", self.dwell["stage1"] * individual.properties("dwell_variation_stage1")

        if individual.properties("trans_prob_stage2") < 0.33:
            yield "crc_pcl2a", self.dwell["stage2"] * individual.properties("dwell_variation_stage2")
            yield "crc_cl2", None
            return
        else:
            yield "crc_pcl2b", self.dwell["stage2"] * individual.properties("dwell_variation_stage2")

        if individual.properties("trans_prob_stage3") < 0.49:
            yield "crc_pcl3a", self.dwell["stage3"] * individual.properties("dwell_variation_stage3")
            yield "crc_cl3", None
            return
        else:
            yield "crc_pcl3b", self.dwell["stage3"] * individual.properties("dwell_variation_stage3")

        yield "crc_pcl4", self.dwell["stage4"] * individual.properties("dwell_variation_stage4")
        yield "crc_cl4", None
        return

    @staticmethod
    def log_durations(individual: Individual, message: str):
        """Method to log the durations."""

        # Only log durations if CRC was clinical
        if "crc_information" not in individual.memory:
            return

        # Obtain the information on the clinical state of CRC of this person
        clinical = individual.memory["crc_information"]

        # Obtain age at which CRC became clinical and age of death
        t_cl = clinical["age"]
        t_death = individual.age

        # Determine clinical stage
        stage = f"{clinical['stage']}"

        # Log time spent in clinical CRC stage
        individual.log_duration(f"crc_cl{stage}", t_cl, t_death)
