from abc import abstractmethod
import numpy as np
from typing import List

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory
from reinvent_scoring.scoring.enums import TransformationTypeEnum


class BasePhysChemComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self._transformation_function = self._assign_transformation(self.parameters.specific_parameters)

    def calculate_score(self, molecules: List) -> ComponentSummary:
        score, raw_score = self._calculate_score(molecules)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters, raw_score=raw_score)
        return score_summary

    def _calculate_score(self, query_mols) -> np.array:
        scores = []
        for mol in query_mols:
            try:
                score = self._calculate_phys_chem_property(mol)
            except ValueError:
                score = 0.0
            scores.append(score)
        transformed_scores = self._transformation_function(scores, self.parameters.specific_parameters)
        return np.array(transformed_scores, dtype=np.float32), np.array(scores, dtype=np.float32)

    @abstractmethod
    def _calculate_phys_chem_property(self, mol):
        raise NotImplementedError("_calculate_phys_chem_property method is not implemented")

    def _assign_transformation(self, specific_parameters: {}):
        transformation_type = TransformationTypeEnum()
        factory = TransformationFactory()
        if self.parameters.specific_parameters[self.component_specific_parameters.TRANSFORMATION]:
            transform_function = factory.get_transformation_function(specific_parameters)
        else:
            self.parameters.specific_parameters[
                self.component_specific_parameters.TRANSFORMATION_TYPE] = transformation_type.NO_TRANSFORMATION
            transform_function = factory.no_transformation
        return transform_function
