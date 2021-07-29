import os
import shutil
import unittest
from unittest.mock import MagicMock

import numpy as np
import numpy.testing as npt

from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.enums import TransformationTypeEnum
from reinvent_scoring.scoring.score_components import RatPKPiP
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH
from unittest_reinvent.scoring_tests.fixtures.predictive_model_fixtures import create_c_lab_component
from unittest_reinvent.scoring_tests.scoring_components.fixtures import score
from unittest_reinvent.fixtures.test_data import CELECOXIB, BUTANE, PENTANE, ETHANOL, ANILINE


class Test_pip_RatPK(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        enum = ScoringFunctionComponentNameEnum()
        csp_enum = ComponentSpecificParametersEnum()
        parameters = create_c_lab_component(enum.RAT_PK_PIP)
        parameters.specific_parameters[csp_enum.TRANSFORMATION] = False
        parameters.specific_parameters[csp_enum.RAT_PK_PROPERTY] = 'log_total_cl_iv'
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, BUTANE, ETHANOL, PENTANE]
        cls.expected_scores = [1.1, 2.4, 1.6, 2.5]

        cls.component = RatPKPiP(parameters)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_pip_1(self):
        result = score(self.component, self.query_smiles)
        has_null = any(np.isnan(result))
        npt.assert_(not has_null, "returned array has Null values")
        unique = np.unique(result)
        npt.assert_equal(np.sort(unique), np.sort(result))

    def test_pip_empty_response(self):
        self.component._parse_single_compound = MagicMock(return_value=np.nan)
        npt.assert_almost_equal(score(self.component, self.query_smiles), [0, 0, 0, 0], 1)


class Test_pip_RatPK_transformation(unittest.TestCase):
    def setUp(cls):
        csp_enum = ComponentSpecificParametersEnum()
        transf_type = TransformationTypeEnum()
        enum = ScoringFunctionComponentNameEnum()
        ts_parameters = create_c_lab_component(enum.RAT_PK_PIP)
        ts_parameters.specific_parameters[csp_enum.RAT_PK_PROPERTY] = 'log_total_cl_iv'
        ts_parameters.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.SIGMOID
        ts_parameters.specific_parameters[csp_enum.HIGH] = 2.4
        ts_parameters.specific_parameters[csp_enum.LOW] = 1
        ts_parameters.specific_parameters[csp_enum.K] = 0.25
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, ANILINE]
        cls.component = RatPKPiP(ts_parameters)

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_pip_transformed_1(self):
        result = score(self.component, self.query_smiles)
        has_null = any(np.isnan(result))
        npt.assert_(not has_null, "returned array has Null values")
        unique = np.unique(result)
        npt.assert_equal(np.sort(unique), np.sort(result))

    def test_pip_empty_response(self):
        self.component._parse_single_compound = MagicMock(return_value=np.nan)
        npt.assert_almost_equal(score(self.component, self.query_smiles), [0, 0], 1)


class Test_pip_RatPK_with_Bioavailability(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        enum = ScoringFunctionComponentNameEnum()
        csp_enum = ComponentSpecificParametersEnum()
        parameters = create_c_lab_component(enum.RAT_PK_PIP)
        parameters.specific_parameters[csp_enum.TRANSFORMATION] = False
        parameters.specific_parameters[csp_enum.RAT_PK_PROPERTY] = 'f_po'
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, BUTANE, ETHANOL, PENTANE]
        cls.expected_scores = [60.1874,22.8223,38.2834, 13.6701]

        cls.component = RatPKPiP(parameters)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_pip_1(self):
        result = score(self.component, self.query_smiles)
        has_null = any(np.isnan(result))
        npt.assert_(not has_null, "returned array has Null values")
        unique = np.unique(result)
        npt.assert_equal(np.sort(unique), np.sort(result))
