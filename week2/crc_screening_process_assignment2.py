from typing import Any, Dict, List

from miscore import Process, Individual


class CRC_screening(Process):
    name = "crc_screening"

    event_tags = [
        "crc_screening_adenoma",  # C
        "crc_screening_cancer",  # C
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

        individual.add_to_queue_age(70, "crc_screening_invitation")

    def crc_screen(self, individual: Individual):
        """Method to process participation and outcomes of invited screening round."""

        # Do nothing if screening was stopped
        if "crc_screening_stop" in individual.memory:
            return

        # B: Implement screening sensitivity
        if individual.memory["crc_lesion"]:
            # B: crc_lesion is the state string itself (e.g. "adenoma"), not a sequence
            stage = individual.memory["crc_lesion"]

            if stage == "adenoma":
                if individual.properties("crc_r_sensitivity_adenoma") < self.sensitivity["adenoma"]:
                    result = "adenoma"  # C: report the lesion type instead of "positive"
                else:
                    result = "negative"

            else:
                if individual.properties("crc_r_sensitivity_all_stage") < self.sensitivity["all_stage"]:
                    result = "cancer"  # C
                else:
                    result = "negative"

            individual.log_event(f"crc_screening_{result}")

        else:
            individual.log_event(f"crc_screening_negative")