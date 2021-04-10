import glob
import os
import re
from abc import abstractmethod
from typing import List

import numpy as np

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory
from reinvent_scoring.scoring.enums import TransformationTypeEnum


class BaseROCSComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self._transformation_function = self._assign_transformation(self.parameters.specific_parameters)

    def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
        return self.calculate_score(molecules, step)

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        # NOTE: valid_idxs are determined with RDKit not with Open Eye
        valid_smiles = self._chemistry.mols_to_smiles(molecules)
        score = self._calculate_omega_score(valid_smiles, step)
        score = self._transformation_function(score, self.parameters.specific_parameters)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters)
        return score_summary

    @abstractmethod
    def _calculate_omega_score(self, smiles, step) -> np.array:
        raise NotImplementedError("_calculate_omega_score method is not implemented")

    def _assign_transformation(self, specific_parameters: {}):
        transformation_type = TransformationTypeEnum()
        factory = TransformationFactory()
        if self.parameters.specific_parameters.get(self.component_specific_parameters.TRANSFORMATION, False):
            transform_function = factory.get_transformation_function(specific_parameters)
        else:
            self.parameters.specific_parameters[
                self.component_specific_parameters.TRANSFORMATION_TYPE] = transformation_type.NO_TRANSFORMATION
            transform_function = factory.no_transformation
        return transform_function

    def _assign_id(self):
        '''
        This is a temp solution until the id=step is passed to this class
        '''
        all_files = glob.glob(os.path.join(self.dir_name, self.overlay_prefix+'*'+'.sdf'))
        if len(all_files) == 0:
            return
        all_files = list(map(os.path.basename, all_files))
        old_ids = list(map(int, [re.findall(r'(\d{4})\.sdf$', f)[0] for f in all_files]))
        self.id = max(old_ids) + 1