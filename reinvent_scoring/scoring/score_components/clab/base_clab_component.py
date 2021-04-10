import json
import subprocess
from abc import abstractmethod
from typing import List

import numpy as np
import rdkit.Chem as rkc

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory
from reinvent_scoring.scoring.enums import TransformationTypeEnum


class BaseClabComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self._clab_input_file = self.parameters.specific_parameters[self.component_specific_parameters.CLAB_INPUT_FILE]
        self.transformation_function = self._assign_transformation(self.parameters.specific_parameters)

    def calculate_score(self, molecules: List[rkc.Mol]) -> ComponentSummary:
        valid_smiles = [rkc.MolToSmiles(mol, isomericSmiles=False) for mol in molecules]
        score, raw_score = self._calculate_score(valid_smiles)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters, raw_score=raw_score)
        return score_summary

    def _calculate_score(self, smiles: List[str]) -> np.array:
        formatted = self._format_clab_input(smiles)
        self._write_clab_input(formatted)
        command = self._create_command(self._clab_input_file)
        results_string = self._send_request_get_stdout(command)
        results_raw = self._parse_multicompound_result(results_string, len(smiles))
        results = self._transformed_results(results_raw)
        results_raw = self._set_nan_to_zero(results_raw)
        return results, results_raw

    def _parse_multicompound_result(self, results: str, data_size: int) -> np.array:
        results_raw = np.empty(data_size, dtype=np.float32)
        results_raw[:] = np.nan
        try:
            loaded = json.loads(results)
            compounds = loaded[0]["Compounds"]
            for compound in compounds:
                try:
                    index = int(compound["Compound"])
                    results_raw[index] = self._parse_compound(compound)
                except (ValueError, TypeError, KeyError):
                    # ValueError: float("a")
                    # TypeError: 3["anything"]
                    # KeyError: compound["bad-name"]
                    pass  # If parsing failed, keep value NaN for this compound and continue.
        finally:
            return results_raw

    def _transformed_results(self, results_raw: np.array) -> np.array:
        """Returns np.array with non-NaN elements transformed by transformation function, and all NaN elements
        transformed into 0. """
        valid_mask = ~np.isnan(results_raw)
        results_raw_valid = results_raw[valid_mask]
        results_transformed = self.transformation_function(results_raw_valid, self.parameters.specific_parameters)
        results = np.zeros(len(results_raw), dtype=np.float32)
        results[valid_mask] = results_transformed
        return results

    def _set_nan_to_zero(self, results_raw: np.array) -> np.array:
        """Returns np.array with  all NaN elements transformed to 0. """
        valid_mask = ~np.isnan(results_raw)
        results_raw_valid = results_raw[valid_mask]
        results = np.zeros(len(results_raw), dtype=np.float32)
        results[valid_mask] = results_raw_valid
        return results

    @abstractmethod
    def _create_command(self, input_file) -> str:
        raise NotImplementedError("_create_command method is not implemented")

    @staticmethod
    def _send_request_get_stdout(command: str) -> str:
        with subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, text=True) as proc:
            stdout, stderr = proc.communicate()
        return stdout

    @abstractmethod
    def _parse_compound(self, compound):
        pass

    def _format_clab_input(self, smiles):
        """Prepares clab input compatible with clab_caco2_efflux, which requires quite strange json.

        Formats input to resemble the following for two molecules:
            {"smiles": "CCC AZ1\nCCCCC AZ2"}

        Alternative input is a simple SMILES (.smi) file, but it would require different extension (.smi).
        Sending simple SMILES-file with .json extension results in an error.
        Corresponding SMILES file example would be:
            CCC AZ1
            CCCCC AZ2
        """

        lines = [f"{smi} {i}" for i, smi in enumerate(smiles)]
        joined_lines = "\\n".join(lines)
        return f'{{"smiles": "{joined_lines}"}}'

    def _write_clab_input(self, data):
        with open(self._clab_input_file, 'w') as ofile:
            ofile.writelines(data)

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
