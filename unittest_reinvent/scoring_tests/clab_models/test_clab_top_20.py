import os
import shutil
import unittest
from unittest.mock import MagicMock

import numpy as np
import numpy.testing as npt

from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH
from unittest_reinvent.fixtures.test_data import CELECOXIB, PHENYLETHANMINE, PROPANE_1_2_DIOL, ACETIC_ACID
from unittest_reinvent.scoring_tests.fixtures.predictive_model_fixtures import create_c_lab_component
from unittest_reinvent.scoring_tests.scoring_components.fixtures import score
from reinvent_scoring.scoring.score_components.clab.top_20 import Top20


class Test_top_20_base(unittest.TestCase):
    def setUp(cls):
        enum = ScoringFunctionComponentNameEnum()
        csp_enum = ComponentSpecificParametersEnum()
        parameters = create_c_lab_component(enum.CLAB_TOP_20)
        parameters.specific_parameters[csp_enum.CLAB_TOP_20_VALUE] = "Base1 pKa"
        parameters.specific_parameters[csp_enum.TRANSFORMATION] = False
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

            cls.query_smiles = [CELECOXIB, PHENYLETHANMINE, PROPANE_1_2_DIOL]
        cls.expected_raw_scores = [0, 10.0, 0]
        cls.component = Top20(parameters)

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_base(self):
        npt.assert_almost_equal(score(self.component, self.query_smiles), self.expected_raw_scores, decimal=1)

    def test_pip_empty_response(self):
        self.component._parse_compound = MagicMock(return_value=np.nan)
        npt.assert_almost_equal(score(self.component, self.query_smiles), [0, 0, 0], 3)


class Test_top_20_acid(unittest.TestCase):
    def setUp(cls):
        enum = ScoringFunctionComponentNameEnum()
        csp_enum = ComponentSpecificParametersEnum()
        parameters = create_c_lab_component(enum.CLAB_TOP_20)
        parameters.specific_parameters[csp_enum.CLAB_TOP_20_VALUE] = "Acid1 pKa"
        parameters.specific_parameters[csp_enum.TRANSFORMATION] = False
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, PHENYLETHANMINE, ACETIC_ACID]
        cls.expected_raw_scores = [9.67, 0.0, 4.79]
        cls.component = Top20(parameters)

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_acid(self):
        npt.assert_almost_equal(score(self.component, self.query_smiles), self.expected_raw_scores, decimal=1)


class Test_top_20_ion_class(unittest.TestCase):
    def setUp(cls):
        enum = ScoringFunctionComponentNameEnum()
        csp_enum = ComponentSpecificParametersEnum()
        parameters = create_c_lab_component(enum.CLAB_TOP_20)
        parameters.specific_parameters[csp_enum.CLAB_TOP_20_VALUE] = "Ion class"
        parameters.specific_parameters["Ion class"] = ['Base', 'Neutral']
        parameters.specific_parameters[csp_enum.TRANSFORMATION] = False
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        cls.query_smiles = [CELECOXIB, PHENYLETHANMINE, ACETIC_ACID]
        cls.expected_raw_scores = [1., 1.0, 0.]
        cls.component = Top20(parameters)

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    def test_ion_class(self):
        npt.assert_almost_equal(score(self.component, self.query_smiles), self.expected_raw_scores, decimal=1)

