from typing import Any, Dict, List

from miscore import Process, Individual


class CRC_screening(Process):
    name = "crc_screening"

    event_tags = [
        "crc_screening_positive",
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

        if individual.memory["crc_lesion"]:
            result = "positive"
        else:
            result = "negative"

        individual.log_event(f"crc_screening_{result}")
