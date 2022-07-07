from joblib import Parallel, delayed
import json
from functools import partial
from typing import List

import numpy as np
import requests
from rdkit import Chem
from requests.adapters import HTTPAdapter
from requests.structures import CaseInsensitiveDict
from requests.packages.urllib3.util.retry import Retry

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary


def batchify(iterable, batch_size):
    for ndx in range(0, len(iterable), batch_size):
        batch = iterable[ndx : min(ndx + batch_size, len(iterable))]
        yield batch


class DlfApiComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.pdb_id = self.parameters.specific_parameters['pdb_id']
        self.chain_id = self.parameters.specific_parameters['chain_id']
        self.dlf_endpoint = self.parameters.specific_parameters.get('dlf_endpoint', "http://ubuntu@3.121.40.10:8080")
        self.wandb_checkpoint = self.parameters.specific_parameters['wandb_checkpoint']
        self.wandb_checkpoint_project = self.parameters.specific_parameters['wandb_checkpoint_project']
        self.batch_size = self.parameters.specific_parameters.get('batch_size', 1)
        self.n_jobs = self.parameters.specific_parameters.get('n_jobs', 4)
        self.session = self.requests_retry_session()

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        score, raw_score = self._get_parallel_score_from_server(molecules)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters, raw_score=raw_score)

        return score_summary

    def _get_parallel_score_from_server(self, query_mols) -> np.array:
        scores = Parallel(n_jobs=self.n_jobs, verbose=10)(
            delayed(self._get_score_from_server)(batch) for batch in batchify(query_mols, batch_size=self.batch_size)
        )
        scores = [item for sublist in scores for item in sublist]
        print(scores)

        transform_params = self.parameters.specific_parameters.get(
            self.component_specific_parameters.TRANSFORMATION, {}
        )
        transformed_scores = self._transformation_function(scores, transform_params)
        return np.array(transformed_scores, dtype=np.float32), np.array(scores, dtype=np.float32)

    def _get_score_from_server(self, query_mols):
        smiles_batch = list(map(partial(Chem.MolToSmiles, canonical=True), query_mols))
        headers = CaseInsensitiveDict()
        headers["Content-Type"] = "application/json"
        data = {
            "smiles": smiles_batch,
            "pdb_id": self.pdb_id,
            "wandb_checkpoint": self.wandb_checkpoint,
            "wandb_checkpoint_project": self.wandb_checkpoint_project,
        }
        json_data = json.dumps(data)
        resp = self.session.post(self.dlf_endpoint, headers=headers, data=json_data)
        if resp.status_code == 200:
            resp = resp.json()
            score_list = resp['output']
        else:
            score_list = []
            print("Retrying one by one...")
            for smiles in smiles_batch:
                data['smiles'] = [smiles]
                json_data = json.dumps(data)
                resp = self.session.post(self.dlf_endpoint, headers=headers, data=json_data)
                if resp.status_code == 200:
                    resp = resp.json()
                    score = resp['output'][0]
                else:
                    print(f"Didn't manage to get a score for {smiles}")
                    score = 0.0
                score_list.append(score)
        return score_list

    @staticmethod
    def requests_retry_session(
            retries=5,
            backoff_factor=1,
            status_forcelist=(500, 502, 504),
            session=None,
    ):
        session = session or requests.Session()
        retry = Retry(
            total=retries,
            read=retries,
            connect=retries,
            backoff_factor=backoff_factor,
            status_forcelist=status_forcelist,
        )
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        return session
