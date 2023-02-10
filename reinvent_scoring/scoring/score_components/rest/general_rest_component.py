import json
import requests

import numpy as np

from typing import List

from reinvent_scoring.scoring.score_components.pip.base_rest_component import BaseRESTComponent
from reinvent_scoring.scoring.component_parameters import ComponentParameters


class GeneralRESTComponent(BaseRESTComponent):
    def __init__(self, parameters: ComponentParameters):
        
        self._server_url = parameters.specific_parameters["server_url"]
        self._server_port = parameters.specific_parameters["server_port"]
        self._server_endpoint = parameters.specific_parameters["server_endpoint"]
        
        self._predictor_id = parameters.specific_parameters["predictor_id"]
        self._predictor_version = parameters.specific_parameters["predictor_version"]
        
        self._request_header = parameters.specific_parameters.get("header", self._default_header)
         
        super().__init__(parameters)
        
    @property
    def _default_header(self):
        return {
                'accept': 'application/json',
                'Content-Type': 'application/json',
            }
        
    def _execute_request(self, request_url, data, header) -> dict:
        params = self._create_params()
        request = requests.post(request_url, json=data, headers=header, params=params)
        if request.status_code != 200:
            raise ValueError(
                f" Status: {request.status_code} Reason: ({request.reason})."
                f"Response content: {request.content}\n"
                f"Response content: {request.text}"
            )
        return request.json()

    def _parse_response(self, response_json: dict, data_size: int) -> np.array:
        compounds = response_json['output']["successes_list"]
        results_raw = np.empty(data_size, dtype=np.float32)
        results_raw[:] = np.nan

        try:
            for compound in compounds:
                try:
                    index = int(compound["query_id"])
                    results_raw[index] = self._parse_single_compound(compound)

                except (ValueError, TypeError, KeyError):
                    pass  # If parsing failed, keep value NaN for this compound and continue.
        finally:
            return results_raw
        
    def _parse_single_compound(self, compound):
        return float(compound["output_value"])
        
    def _format_data(self, smiles: List[str]) -> dict:
        json_data = [{'input_string': smi,
                      'query_id': str(i)} for i, smi in enumerate(smiles)]
        return json_data
    
    def _create_url(self, component_name) -> str:
        url = f"{self._server_url}:{self._server_port}/{self._server_endpoint}"
        return url
    
    
    def _create_header(self) -> dict:
        return self._request_header
    
    def _create_params(self) -> dict:
        return {
                'predictor_id': self._predictor_id,
                'predictor_version': self._predictor_version,
                'inp_fmt': 'smiles',
            }