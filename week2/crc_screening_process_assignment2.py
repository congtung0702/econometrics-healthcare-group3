from typing import Any, Dict, List

from miscore import Process, Individual
from crc_process_assignment2 import CANCER_STAGE  # E


class CRC_screening(Process):
    name = "crc_screening"

    # C: Change positive result into adenoma or cancer stage tests
    event_tags = [
        "crc_screening_adenoma",
        "crc_screening_cancer",
        "crc_screening_negative",
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
    def stop_screening(individual: Individual):
        individual.memory["crc_screening_stop"] = True

    def schedule_screening(self, individual: Individual):
        """Method to schedule the screening invitations."""

        # H: Add screening ages at correct ages
        for age in self.screening_ages:
            individual.add_to_queue_age(age, "crc_screening_invitation")

    def crc_screen(self, individual: Individual):
        """Method to process participation and outcomes of invited screening round."""

        # Do nothing if screening was stopped
        if "crc_screening_stop" in individual.memory:
            return

        # B: Implement screening sensitivity
        if individual.memory["crc_lesion"]:
            stage = individual.memory["crc_lesion"]

            # D: screen-detected adenomas become known, and are removed immediately
            if stage == "adenoma":
                if individual.properties("crc_r_sensitivity_adenoma") < self.sensitivity["adenoma"]:
                    # C: return adenoma stage
                    result = "adenoma"
                    individual.memory["crc_lesion_removed"] = True
                    individual.memory["crc_lesion"] = []
                else:
                    result = "negative"

            # E: screen-detected cancer becomes clinically known, at the corresponding stage
            else:
                if individual.properties("crc_r_sensitivity_all_stage") < self.sensitivity["all_stage"]:
                    # C: return cancer stage
                    result = "cancer"
                    # F: Pass screen_detected boolean as TRUE
                    individual.add_to_queue(0., "crc_clinical", stage=CANCER_STAGE[stage], screen_detected=True)
                else:
                    result = "negative"

            individual.log_event(f"crc_screening_{result}")

        else:
            individual.log_event(f"crc_screening_negative")