import json
import os
from abc import abstractmethod
from typing import List

import numpy as np

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.enums import EnvironmentalVariablesEnum
from reinvent_scoring.scoring.score_components.pip.base_rest_component import BaseRESTComponent


class BasePiPModelComponent(BaseRESTComponent):
    def __init__(self, parameters: ComponentParameters):
        self._environment_keys = EnvironmentalVariablesEnum()
        super().__init__(parameters)

    def _parse_response(self, response_json: dict, data_size: int) -> np.array:
        compounds = response_json['jsonData']['data']
        results_raw = np.empty(data_size, dtype=np.float32)
        results_raw[:] = np.nan

        try:
            for compound in compounds:
                try:
                    index = int(compound["id"])
                    results_raw[index] = self._parse_single_compound(compound)

                except (ValueError, TypeError, KeyError):
                    pass  # If parsing failed, keep value NaN for this compound and continue.
        finally:
            return results_raw

    @abstractmethod
    def _parse_single_compound(self, compound):
        raise NotImplementedError("_parse_compound method is not implemented")

    def _format_data(self, smiles: List[str]) -> dict:
        molecules = [{"molData": smi, "id": f"{i}"} for i, smi in enumerate(smiles)]
        data = {
            "jsonData": {
                "data": molecules,
                "metadata": {
                    "molFormat":
                        "smiles"
                },
                "parameters": {}
            }
        }
        return data

    def _create_url(self, component_name) -> str:
        pip_url = self._get_enviornment_variable(self._environment_keys.PIP_URL)
        request_url = pip_url.format(component_name)
        return request_url

    def _create_header(self) -> dict:
        pip_key = self._get_enviornment_variable(self._environment_keys.PIP_KEY)

        header = {'Content-Type': 'application/json', 'x-api-key': pip_key}
        return header

    def _get_enviornment_variable(self, variable: str) -> str:
        try:
            return os.environ[variable]
        except KeyError:
            return self._retrieve_pip_key_from_config(variable)

    def _retrieve_pip_key_from_config(self, variable: str) -> str:
        try:
            project_root = os.path.dirname(__file__)
            with open(os.path.join(project_root, '../../../configs/config.json'), 'r') as f:
                config = json.load(f)
            environmental_variables = config[self._environment_keys.ENVIRONMENTAL_VARIABLES]
            return environmental_variables[variable]
        except KeyError as ex:
            raise KeyError(f"Key {variable} not found in reinvent_scoring/configs/config.json")