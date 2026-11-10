from typing import Any, Dict, List
import random

from miscore import Process, Individual

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
        "crc_screening_positive",
        "crc_screening_negative",
        "crc_screening_adenoma",
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

    @staticmethod
    def stop_screening(Individual: Individual):
        Individual.memory["crc_screening_stop"] = True

    def schedule_screening(self, Individual: Individual):
        """Method to schedule the screening invitations."""
        for scr_age in self.screening_ages:
            Individual.add_to_queue_age(scr_age, "crc_screening_invitation")

    def crc_screen(self, Individual: Individual):
        """Method to process participation and outcomes of invited screening round."""

        # Do nothing if screening was stopped
        if "crc_screening_stop" in Individual.memory:
            return

        if Individual.memory["crc_lesion"]:
            if random.random() < self.sensitivity[Individual.memory["crc_lesion"]]:
                if Individual.memory["crc_lesion"] == "adenoma":
                    result = "adenoma"
                    Individual.memory["crc_lesion_removed"] = True
                else:
                    result = "cancer"
                    Individual.add_to_queue(0., "crc_clinical",
                                         stage=CANCER_STAGE[Individual.memory["crc_lesion"]],
                                         screen_detected=True)

            else:
                result = "negative"
        else:
            result = "negative"

        Individual.log_event(f"crc_screening_{result}")
