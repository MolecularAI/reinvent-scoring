import json
import os
import shutil
import unittest

import numpy.testing as npt
import pytest

from reinvent_chemistry.conversions import Conversions
from reinvent_scoring.scoring.score_components.console_invoked.azgard import AZgard
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH, AZGARD_UNITTEST_JSON, AZGARD_UNITTEST_NIBR_NEGATIVE_IMAGE, \
                                             AZGARD_UNITTEST_GRID_PATH, AZGARD_UNITTEST_NIBR_VALUES_KEY
from unittest_reinvent.scoring_tests.fixtures.predictive_model_fixtures import create_AZgard_component_parameters
from unittest_reinvent.fixtures.test_data import PARACETAMOL, ASPIRIN, CAFFEINE
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from unittest_reinvent.unit_testing import ignore_warnings


class Test_structural_AZgard(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.csp_enum = ComponentSpecificParametersEnum()
        chemistry = Conversions()
        if not os.path.isdir(MAIN_TEST_PATH):
            os.makedirs(MAIN_TEST_PATH)
        cls.old_dir = os.getcwd()
        os.chdir(MAIN_TEST_PATH)

        # load internal "AZgard" configuration file and update the receptor and negative-image file paths, respectively
        with open(AZGARD_UNITTEST_JSON, 'r') as f:
            azgard_conf = json.load(f)
        azgard_conf["workflow"]["steps"][2]["settings"]["additional"]["configuration"]["GRIDFILE"] = [AZGARD_UNITTEST_GRID_PATH]
        azgard_conf["workflow"]["steps"][3]["input"]["generic"][0]["source"] = AZGARD_UNITTEST_NIBR_NEGATIVE_IMAGE

        # write a temporary copy of the configuration file out
        temp_conf_path = MAIN_TEST_PATH + "/temp_azgard_NIBR_unittest.json"
        with open(temp_conf_path, 'w') as f:
            json.dump(azgard_conf, f, indent=4)

        # initialize AZgard with three random molecules and update path to temporary file
        parameters = create_AZgard_component_parameters()
        parameters.specific_parameters[cls.csp_enum.AZGARD_CONFPATH] = temp_conf_path
        parameters.specific_parameters[cls.csp_enum.AZGARD_VALUES_KEY] = AZGARD_UNITTEST_NIBR_VALUES_KEY
        cls.component = AZgard(parameters)
        cls.query_smiles = [ASPIRIN, PARACETAMOL, CAFFEINE]
        cls.query_mols = [chemistry.smile_to_mol(smile) for smile in cls.query_smiles]

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.old_dir)
        if os.path.isdir(MAIN_TEST_PATH):
            shutil.rmtree(MAIN_TEST_PATH)

    @pytest.mark.integration
    @ignore_warnings
    def test_AZgard(self):
        summary = self.component.calculate_score(self.query_mols)
        npt.assert_almost_equal(summary.total_score, [0.476677, 0.458017, 0.510676], decimal=3)
        self.assertEqual(len(summary.total_score), 3)
