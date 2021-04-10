import io
import subprocess
from abc import abstractmethod

import numpy as np
from typing import List
import rdkit.Chem as rkc

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory
from reinvent_scoring.scoring.enums import TransformationTypeEnum


class BaseStructuralComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.transformation_function = self._assign_transformation(self.parameters.specific_parameters)

    def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
        return self.calculate_score(molecules, step)

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        valid_smiles = [rkc.MolToSmiles(mol, isomericSmiles=False) for mol in molecules]
        score, raw_score = self._calculate_score(valid_smiles, step)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters, raw_score=raw_score)
        return score_summary

    @abstractmethod
    def _calculate_score(self, smiles: List[str], step) -> np.array:
        raise NotImplementedError("_calculate_score method is not implemented")

    @abstractmethod
    def _create_command(self, input_file, step) -> str:
        raise NotImplementedError("_create_command method is not implemented")

    def _send_request_with_stepwize_read(self, command, data_size: int):
        with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              shell=True) as proc:
            wrapped_proc_in = io.TextIOWrapper(proc.stdin, 'utf-8')
            wrapped_proc_out = io.TextIOWrapper(proc.stdout, 'utf-8')
            result = [self._parse_result(wrapped_proc_out.readline()) for i in range(data_size)]
            wrapped_proc_in.close()
            wrapped_proc_out.close()
            proc.wait()
            proc.terminate()
        return result

    @abstractmethod
    def _parse_result(self, result):
        raise NotImplementedError("_parse_result method is not implemented")

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

