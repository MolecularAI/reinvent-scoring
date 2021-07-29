import unittest

from reinvent_scoring import ScoringFunctionComponentNameEnum, ComponentParameters


class BaseTestMatchingSubstructure(unittest.TestCase):

    def setUp(self):
        sf_enum = ScoringFunctionComponentNameEnum()
        self.parameters = ComponentParameters(component_type=sf_enum.MATCHING_SUBSTRUCTURE,
                                              name="matching_substructure",
                                              weight=1.,
                                              smiles=self.smiles,
                                              model_path="",
                                              specific_parameters={})
