import unittest
from typing import List

from reinvent_scoring import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring import CustomSum


class BaseTestCustomSum(unittest.TestCase):

    def setup_attrs(self):
        self.sf_enum = ScoringFunctionComponentNameEnum()

    def init(self, component_type: str, name: str, smiles: List[str] = []):
        self.component_type = component_type
        self.name = name
        self.smiles = smiles

    def setUp(self):
        self.parameters = ComponentParameters(component_type=self.component_type,
                                              name=self.name,
                                              weight=1.,
                                              smiles=self.smiles,
                                              model_path="",
                                              specific_parameters={})
        self.sf_state = CustomSum(parameters=[self.parameters])
