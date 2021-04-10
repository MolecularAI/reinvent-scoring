import asyncio
from abc import abstractmethod
from typing import List

import aiohttp
# import nest_asyncio
import numpy as np
import rdkit.Chem as rkc

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.enums import TransformationTypeEnum
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory


# nest_asyncio.apply()


class BaseRESTComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        # self._environment_keys = EnvironmentalVariablesEnum()
        self.transformation_function = self._assign_transformation(self.parameters.specific_parameters)
        self._request_url = self._create_url(self.parameters.component_type)
        self._request_header = self._create_header()

    def calculate_score(self, molecules: List[rkc.Mol]) -> ComponentSummary:
        valid_smiles = self._chemistry.mols_to_smiles(molecules)
        score, raw_score = self._score_smiles(valid_smiles)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters, raw_score=raw_score)

        return score_summary

    def _score_smiles(self, smiles: List[str]) -> np.array:
        loop = asyncio.get_event_loop()
        response = loop.run_until_complete(self._post_request(self._request_url, smiles, self._request_header))
        results_raw = self._parse_response(response, len(smiles))
        results = self._apply_score_transformation(results_raw)

        return results, results_raw

    async def _post_request(self, url, smiles, header):
        data = self._format_data(smiles)
        result = await self._execute_request(url, data, header)

        return result

    async def _execute_request(self, request_url, data, header) -> str:
        async with aiohttp.ClientSession() as session:
            async with session.post(request_url, json=data, headers=header) as response:
                if response.status != 200:
                    raise ValueError(
                        f" Status: {response.status} Reason: ({response.reason})."
                        f"Response content: {response.content}"
                    )
                return await response.json()

    @abstractmethod
    def _parse_response(self, response_json: dict, data_size: int) -> np.array:
        raise NotImplementedError("_parse_response method is not implemented")

    def _apply_score_transformation(self, results_raw: np.array) -> np.array:
        """Returns np.array with non-NaN elements transformed by transformation function, and all NaN elements
        transformed into 0. """
        valid_mask = ~np.isnan(results_raw)
        results_raw_valid = results_raw[valid_mask]
        results_transformed = self.transformation_function(results_raw_valid, self.parameters.specific_parameters)
        results = np.zeros(len(results_raw), dtype=np.float32)
        results[valid_mask] = results_transformed

        return results

    @abstractmethod
    def _parse_single_compound(self, compound):
        raise NotImplementedError("_parse_compound method is not implemented")

    @abstractmethod
    def _format_data(self, smiles: List[str]) -> dict:
        raise NotImplementedError("_format_data method is not implemented")

    @abstractmethod
    def _create_url(self, component_name) -> str:
        raise NotImplementedError("_create_url method is not implemented")

    @abstractmethod
    def _create_header(self) -> dict:
        raise NotImplementedError("_create_header method is not implemented")

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
