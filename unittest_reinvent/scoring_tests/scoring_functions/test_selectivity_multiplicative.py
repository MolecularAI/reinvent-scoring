import unittest

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring import CustomProduct
from unittest_reinvent.scoring_tests.fixtures.predictive_model_fixtures import create_activity_component_regression, \
    create_predictive_property_component_regression, create_custom_alerts_configuration
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from unittest_reinvent.fixtures.test_data import CELECOXIB, ASPIRIN,BENZENE, ANILINE, METAMIZOLE, GENTAMICIN, PROPANE


class Test_selectivity_multiplicative_function(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        enum = ScoringFunctionComponentNameEnum()
        activity = create_activity_component_regression()
        qed_score = ComponentParameters(component_type=enum.QED_SCORE,
                                        name="qed_score_name",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})

        custom_alerts = create_custom_alerts_configuration()

        matching_substructure = ComponentParameters(component_type=enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure_name",
                                                    weight=1.,
                                                    smiles=[BENZENE],
                                                    model_path="",
                                                    specific_parameters={})
        self.sf_state = CustomProduct(
            parameters=[activity, qed_score, custom_alerts, matching_substructure])

    def test_special_selectivity_multiplicative_1(self):
        score = self.sf_state.get_final_score(smiles=[ANILINE])
        self.assertAlmostEqual(score.total_score[0], 0.357, 3)

    def test_special_selectivity_multiplicative_2(self):
        score_1 = self.sf_state.get_final_score(smiles=[CELECOXIB])
        score_2 = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score_1.total_score[0], 0.481, 3)
        self.assertAlmostEqual(score_2.total_score[0], 0.481, 3)


class Test_selectivity_multiplicative_function_with_predictive_property(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        enum = ScoringFunctionComponentNameEnum()
        predictive_property = create_predictive_property_component_regression()
        activity = create_activity_component_regression()
        qed_score = ComponentParameters(component_type=enum.QED_SCORE,
                                        name="qed_score_name",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})

        custom_alerts = create_custom_alerts_configuration()

        matching_substructure = ComponentParameters(component_type=enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure_name",
                                                    weight=1.,
                                                    smiles=[BENZENE],
                                                    model_path="",
                                                    specific_parameters={})
        self.sf_state = CustomProduct(
            parameters=[activity, qed_score, custom_alerts, matching_substructure, predictive_property])

    def test_special_selectivity_multiplicative_1(self):
        score = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score.total_score[0], 0.414, 3)


class Test_selectivity_multiplicative_function_with_predictive_property_and_alert(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        enum = ScoringFunctionComponentNameEnum()
        predictive_property = create_predictive_property_component_regression()
        self.activity = create_activity_component_regression()
        qed_score = ComponentParameters(component_type=enum.QED_SCORE,
                                        name="qed_score_name",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})

        custom_alerts = create_custom_alerts_configuration()

        matching_substructure = ComponentParameters(component_type=enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure_name",
                                                    weight=1.,
                                                    smiles=[BENZENE],
                                                    model_path="",
                                                    specific_parameters={})
        self.sf_state = CustomProduct(
            parameters=[self.activity, qed_score, custom_alerts, matching_substructure, predictive_property])

    def test_special_selectivity_multiplicative_with_alert_1(self):
        score = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score.total_score[0], 0.414, 3)

    def test_special_selectivity_multiplicative_with_alert_2(self):
        score = self.sf_state.get_final_score(smiles=[METAMIZOLE])
        self.assertEqual(score.total_score, 0)

    def test_special_selectivity_multiplicative_with_alert_3(self):
        score = self.sf_state.get_final_score(smiles=[ANILINE])
        self.assertAlmostEqual(score.total_score[0], 0.324, 3)

    def test_special_selectivity_multiplicative_with_alert_4(self):
        score = self.sf_state.get_final_score(smiles=[GENTAMICIN])
        self.assertEqual(score.total_score[0], 0)


class Test_selectivity_multiplicative_function_with_predictive_property_no_sigm_trans(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        enum = ScoringFunctionComponentNameEnum()
        predictive_property = create_predictive_property_component_regression()
        activity = create_activity_component_regression()
        qed_score = ComponentParameters(component_type=enum.QED_SCORE,
                                        name="qed_score_name",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})
        custom_alerts = ComponentParameters(component_type=enum.CUSTOM_ALERTS,
                                            name="custom_alerts_name",
                                            weight=1.,
                                            smiles=[PROPANE],
                                            model_path="",
                                            specific_parameters={})
        matching_substructure = ComponentParameters(component_type=enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure_name",
                                                    weight=1.,
                                                    smiles=[BENZENE],
                                                    model_path="",
                                                    specific_parameters={})
        self.sf_state = CustomProduct(
            parameters=[activity, qed_score, custom_alerts, matching_substructure, predictive_property])

    def test_special_selectivity_multiplicative_no_sigm_trans_1(self):
        score = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score.total_score[0], 0.414, 3)

    def test_special_selectivity_multiplicative_no_sigm_trans_2(self):
        score = self.sf_state.get_final_score(smiles=[ASPIRIN])
        self.assertAlmostEqual(score.total_score[0], 0.412, 3)

    def test_special_selectivity_multiplicative_no_sigm_trans_4(self):
        score = self.sf_state.get_final_score(smiles=[GENTAMICIN])
        self.assertEqual(score.total_score[0], 0)
