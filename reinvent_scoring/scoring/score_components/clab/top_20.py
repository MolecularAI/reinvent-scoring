from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.clab.base_clab_component import BaseClabComponent


class Top20(BaseClabComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

    def _create_command(self, input_file) -> str:
        name = "clab_top20"
        command = f"ml {name}; {name} {input_file}"
        return command

    def _parse_compound(self, compound) -> float:
        top_20_value = self.parameters.specific_parameters[self.component_specific_parameters.CLAB_TOP_20_VALUE]
        if top_20_value == self.component_specific_parameters.ION_CLASS:
            ion_class = self.parameters.specific_parameters[top_20_value]
            result = any(compound[top_20_value] == c for c in ion_class)
            return float(result)
        return float(compound[top_20_value])
