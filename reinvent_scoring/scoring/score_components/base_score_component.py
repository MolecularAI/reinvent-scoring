from abc import ABC, abstractmethod
from typing import List

from reinvent_scoring.chemistry.general import GeneralChemistry
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum


class BaseScoreComponent(ABC):

    def __init__(self, parameters: ComponentParameters):
        self.component_specific_parameters = ComponentSpecificParametersEnum()
        self.parameters = parameters
        self._chemistry = GeneralChemistry()

    @abstractmethod
    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        raise NotImplementedError("calculate_score method is not implemented")

    def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
        return self.calculate_score(molecules)
