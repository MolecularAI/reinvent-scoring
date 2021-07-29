import os
import shutil
import unittest

from reinvent_scoring.scoring.scoring_function_factory import ScoringFunctionFactory
from reinvent_scoring.scoring.scoring_function_parameters import ScoringFunctionParameters
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH
from unittest_reinvent.fixtures.test_data import CELECOXIB
from unittest_reinvent.scoring_tests.fixtures import create_c_lab_component, \
    create_predictive_property_component_regression
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.enums import ScoringFunctionNameEnum
from unittest_reinvent.unit_testing import ignore_warnings


class Test_pip_with_parallel_additive(unittest.TestCase):
    def setUp(self):
        csp_enum = ComponentSpecificParametersEnum()
        sf_enum = ScoringFunctionComponentNameEnum()
        sf_name_enum = ScoringFunctionNameEnum()

        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)

        activity = create_predictive_property_component_regression()

        solph74 = create_c_lab_component(sf_enum.SOLUBILITY_DD_PIP)
        solph74.weight = 4
        solph74.specific_parameters[csp_enum.HIGH] = 3
        solph74.specific_parameters[csp_enum.LOW] = 1

        hlm_clint = create_c_lab_component(sf_enum.HLM_CLINT_PIP)
        hlm_clint.weight = 4
        hlm_clint.specific_parameters[csp_enum.HIGH] = 3
        hlm_clint.specific_parameters[csp_enum.LOW] = 1

        azlogd74 = create_c_lab_component(sf_enum.AZ_LOGD74_PIP)
        azlogd74.weight = 4
        azlogd74.specific_parameters[csp_enum.HIGH] = 3
        azlogd74.specific_parameters[csp_enum.LOW] = 1

        sf_parameters = ScoringFunctionParameters(name=sf_name_enum.CUSTOM_SUM,
                                                  parameters=[vars(hlm_clint), vars(activity), vars(solph74),
                                                              vars(azlogd74)], parallel=True)
        self.sf_instance = ScoringFunctionFactory(sf_parameters=sf_parameters)

    def tearDown(self):
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    @ignore_warnings
    def test_parallel_clab_1(self):
        smiles = [CELECOXIB]*3
        score = self.sf_instance.get_final_score(smiles=smiles)
        self.assertAlmostEqual(score.total_score[0], 0.348, 2)
