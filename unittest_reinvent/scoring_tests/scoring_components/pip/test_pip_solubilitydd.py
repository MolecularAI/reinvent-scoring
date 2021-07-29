import os
import shutil
import unittest
from unittest.mock import MagicMock

import numpy as np
import numpy.testing as npt

from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.enums import TransformationTypeEnum
from reinvent_scoring.scoring.score_components.pip.pip_log_prediction_component import PiPLogPredictionComponent
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH
from unittest_reinvent.scoring_tests.fixtures.predictive_model_fixtures import create_c_lab_component
from unittest_reinvent.scoring_tests.scoring_components.fixtures import score
from unittest_reinvent.fixtures.test_data import CELECOXIB, BUTANE, ETHANOL, BUTAN_1_AMINE


class Test_pip_solubilityDD(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        enum = ScoringFunctionComponentNameEnum()
        csp_enum = ComponentSpecificParametersEnum()
        parameters = create_c_lab_component(enum.SOLUBILITY_DD_PIP)
        parameters.specific_parameters[csp_enum.TRANSFORMATION] = False
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, BUTANE, ETHANOL, BUTAN_1_AMINE]
        cls.expected_scores = [0.282, 1.071, 1.713,  1.628]

        cls.component = PiPLogPredictionComponent(parameters)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_pip_1(self):
        result = score(self.component, self.query_smiles)
        # result[0] = np.nan
        has_null = any(np.isnan(result))
        npt.assert_(not has_null, "returned array has Null values")
        unique = np.unique(result)
        npt.assert_equal(np.sort(unique), np.sort(result))

    def test_pip_empty_response(self):
        self.component._parse_single_compound = MagicMock(return_value=np.nan)
        npt.assert_almost_equal(score(self.component, self.query_smiles), [0, 0, 0, 0], 3)


class Test_pip_solubilityDD_transformation(unittest.TestCase):
    def setUp(cls):
        csp_enum = ComponentSpecificParametersEnum()
        transf_type = TransformationTypeEnum()
        enum = ScoringFunctionComponentNameEnum()
        ts_parameters = create_c_lab_component(enum.SOLUBILITY_DD_PIP)
        ts_parameters.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.STEP
        ts_parameters.specific_parameters[csp_enum.HIGH] = 3
        ts_parameters.specific_parameters[csp_enum.LOW] = 1
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, BUTANE, ETHANOL, BUTAN_1_AMINE]
        cls.expected_scores = [0.0, 0.0, 1.0, 1.0]
        cls.component = PiPLogPredictionComponent(ts_parameters)

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_pip_transformed_1(self):
        npt.assert_almost_equal(score(self.component, self.query_smiles), self.expected_scores, decimal=3)

    def test_pip_empty_response(self):
        self.component._parse_single_compound = MagicMock(return_value=np.nan)
        npt.assert_almost_equal(score(self.component, self.query_smiles), [0, 0, 0, 0], 3)
