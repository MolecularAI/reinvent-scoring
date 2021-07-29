import os
import shutil
import unittest

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring import CustomProduct
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH
from unittest_reinvent.scoring_tests.fixtures import create_c_lab_component, \
    create_predictive_property_component_regression
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.enums import TransformationTypeEnum
from unittest_reinvent.fixtures.test_data import CELECOXIB


class Test_pip_with_Product(unittest.TestCase):
    def setUp(self):
        csp_enum = ComponentSpecificParametersEnum()
        sf_enum = ScoringFunctionComponentNameEnum()
        transf_type = TransformationTypeEnum()

        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        qed_score = ComponentParameters(component_type=sf_enum.QED_SCORE,
                                        name="qed_score_name",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})
        custom_alerts = ComponentParameters(component_type=sf_enum.CUSTOM_ALERTS,
                                            name="custom_alerts_name",
                                            weight=1.,
                                            smiles=[],
                                            model_path="",
                                            specific_parameters={})
        matching_substructure = ComponentParameters(component_type=sf_enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure_name",
                                                    weight=1.,
                                                    smiles=[],
                                                    model_path="",
                                                    specific_parameters={})

        solDD = create_c_lab_component(sf_enum.SOLUBILITY_DD_PIP)
        solDD.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.SIGMOID
        solDD.specific_parameters[csp_enum.K] = 0.5
        solDD.specific_parameters[csp_enum.LOW] = 0.3
        solDD.specific_parameters[csp_enum.HIGH] = 2.3
        solDD.specific_parameters[csp_enum.TRANSFORMATION] = True
        solDD.weight = 4

        activity = create_predictive_property_component_regression()

        self.sf_state = CustomProduct(parameters=[solDD, qed_score, custom_alerts, matching_substructure, activity])

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_clab_with_primary_multiplicative_1(self):
        smiles = [CELECOXIB]*3
        score = self.sf_state.get_final_score(smiles=smiles)
        self.assertAlmostEqual(score.total_score[0], 0.025, 3)
