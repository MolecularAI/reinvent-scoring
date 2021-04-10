from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.clab.base_clab_component import BaseClabComponent


class RHClint(BaseClabComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

    def _create_command(self, input_file):
        name = "clab_rh_clint"
        command = f"ml {name}; {name} {input_file}"
        return command

    def _parse_compound(self, compound):
        return float(compound["Prediction"]["point_estimate"])
