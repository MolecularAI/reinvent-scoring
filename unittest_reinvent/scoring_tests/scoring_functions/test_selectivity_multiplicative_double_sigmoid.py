import unittest

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring import CustomProduct
from reinvent_scoring.scoring.score_summary import FinalSummary
from unittest_reinvent.fixtures.paths import ACTIVITY_REGRESSION
from unittest_reinvent.scoring_tests.fixtures.predictive_model_fixtures import create_activity_component_regression, \
    create_offtarget_activity_component_regression, create_custom_alerts_configuration
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.enums import TransformationTypeEnum
from unittest_reinvent.fixtures.test_data import CELECOXIB, ASPIRIN, AMOXAPINE, HEXANE, BUTANE, ANILINE, GENTAMICIN, \
    BENZENE


class Test_desirability_multiplicative_function(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        csp_enum = ComponentSpecificParametersEnum()
        transf_type = TransformationTypeEnum()
        sf_enum = ScoringFunctionComponentNameEnum()
        activity = create_activity_component_regression()
        activity.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.DOUBLE_SIGMOID
        activity.specific_parameters[csp_enum.COEF_DIV] = 100.
        activity.specific_parameters[csp_enum.COEF_SI] = 150.
        activity.specific_parameters[csp_enum.COEF_SE] = 150.
        off_activity = create_offtarget_activity_component_regression()

        delta_params = {
            "high": 3.0,
            "k": 0.25,
            "low": 0.0,
            "transformation": True,
            "transformation_type": "sigmoid"
        }

        selectivity = ComponentParameters(component_type=sf_enum.SELECTIVITY,
                                           name="desirability",
                                           weight=1.,
                                           smiles=[],
                                           model_path="",
                                           specific_parameters={
                                               "activity_model_path": activity.model_path,
                                               "offtarget_model_path": off_activity.model_path,
                                               "activity_specific_parameters": activity.specific_parameters.copy(),
                                               "offtarget_specific_parameters": off_activity.specific_parameters,
                                               "delta_transformation_parameters": delta_params
                                           })

        qed_score = ComponentParameters(component_type=sf_enum.QED_SCORE,
                                        name="qed_score",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})
        matching_substructure = ComponentParameters(component_type=sf_enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure",
                                                    weight=1.,
                                                    smiles=[BENZENE],
                                                    model_path="",
                                                    specific_parameters={})

        custom_alerts = create_custom_alerts_configuration()

        self.sf_state = CustomProduct(
            parameters=[activity, selectivity, qed_score, matching_substructure, custom_alerts])

    def test_desirability_multiplicative_1(self):
        score = self.sf_state.get_final_score(smiles=[BUTANE])
        self.assertAlmostEqual(score.total_score[0], 0.142, 3)

    def test_desirability_multiplicative_2(self):
        score = self.sf_state.get_final_score(smiles=[ANILINE])
        self.assertAlmostEqual(score.total_score[0], 0.294, 3)

    def test_desirability_multiplicative_3(self):
        score = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score.total_score[0], 0.342, 3)

    def test_desirability_multiplicative_4(self):
        score = self.sf_state.get_final_score(smiles=[HEXANE, "12"])
        self.assertAlmostEqual(score.total_score[0], 0.145, 3)
        self.assertEqual(score.total_score[1], 0)


class Test_special_desirability_multiplicative_function(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        csp_enum = ComponentSpecificParametersEnum()
        transf_type = TransformationTypeEnum()
        sf_enum = ScoringFunctionComponentNameEnum()
        activity = create_activity_component_regression()
        activity.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.DOUBLE_SIGMOID
        activity.specific_parameters[csp_enum.COEF_DIV] = 100.
        activity.specific_parameters[csp_enum.COEF_SI] = 150.
        activity.specific_parameters[csp_enum.COEF_SE] = 150.
        off_activity = create_offtarget_activity_component_regression()
        off_activity.weight = 4
        off_activity.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.DOUBLE_SIGMOID
        off_activity.specific_parameters[csp_enum.COEF_DIV] = 100.
        off_activity.specific_parameters[csp_enum.COEF_SI] = 150.
        off_activity.specific_parameters[csp_enum.COEF_SE] = 150.
        off_activity.model_path = ACTIVITY_REGRESSION
        off_activity2 = create_offtarget_activity_component_regression()
        off_activity2.weight = 4

        delta_params = {
            "high": 3.0,
            "k": 0.25,
            "low": 0.0,
            "transformation": True,
            "transformation_type": "sigmoid"
        }

        selectivity_1 = ComponentParameters(component_type=sf_enum.SELECTIVITY,
                                           name="selectivity_1",
                                           weight=4.,
                                           smiles=[],
                                           model_path="",
                                           specific_parameters={
                                               "activity_model_path": activity.model_path,
                                               "offtarget_model_path": off_activity.model_path,
                                               "activity_specific_parameters": activity.specific_parameters.copy(),
                                               "offtarget_specific_parameters": off_activity.specific_parameters,
                                               "delta_transformation_parameters": delta_params
                                           })

        selectivity_2 = ComponentParameters(component_type=sf_enum.SELECTIVITY,
                                           name="selectivity_2",
                                           weight=4.,
                                           smiles=[],
                                           model_path="",
                                           specific_parameters={
                                               "activity_model_path": activity.model_path,
                                               "offtarget_model_path": off_activity2.model_path,
                                               "activity_specific_parameters": activity.specific_parameters.copy(),
                                               "offtarget_specific_parameters": off_activity2.specific_parameters,
                                               "delta_transformation_parameters": delta_params
                                           })

        qed_score = ComponentParameters(component_type=sf_enum.QED_SCORE,
                                        name="qed_score",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})
        custom_alerts = ComponentParameters(component_type=sf_enum.CUSTOM_ALERTS,
                                            name="custom_alerts",
                                            weight=1.,
                                            smiles=[],
                                            model_path="",
                                            specific_parameters={})
        matching_substructure = ComponentParameters(component_type=sf_enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure",
                                                    weight=1.,
                                                    smiles=[GENTAMICIN],
                                                    model_path="",
                                                    specific_parameters={})
        self.sf_state = CustomProduct(
            parameters=[activity, selectivity_1, selectivity_2, qed_score, custom_alerts, matching_substructure])

    def test_special_desirability_multiplicative_1(self):
        score = self.sf_state.get_final_score(smiles=[ASPIRIN])
        self.assertAlmostEqual(score.total_score[0], 0.0451, 3)

    def test_special_desirability_multiplicative_2(self):
        score = self.sf_state.get_final_score(smiles=[ANILINE])
        self.assertAlmostEqual(score.total_score[0], 0.044, 3)

    def test_special_desirability_multiplicative_3(self):
        score = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score.total_score[0], 0.047, 3)


class Test_selectivity_function_with_double_sigmoid(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        csp_enum = ComponentSpecificParametersEnum()
        transf_type = TransformationTypeEnum()
        enum = ScoringFunctionComponentNameEnum()
        activity = create_activity_component_regression()
        activity.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.DOUBLE_SIGMOID
        activity.specific_parameters[csp_enum.COEF_DIV] = 100.
        activity.specific_parameters[csp_enum.COEF_SI] = 150.
        activity.specific_parameters[csp_enum.COEF_SE] = 150.
        off_activity = create_offtarget_activity_component_regression()

        delta_params = {
            "high": 3.0,
            "k": 0.25,
            "low": 0.0,
            "transformation": True,
            "transformation_type": "sigmoid"
        }

        selectivity = ComponentParameters(component_type=enum.SELECTIVITY,
                                           name="desirability",
                                           weight=1.,
                                           smiles=[],
                                           model_path="",
                                           specific_parameters={
                                               "activity_model_path": activity.model_path,
                                               "offtarget_model_path": off_activity.model_path,
                                               "activity_specific_parameters": activity.specific_parameters.copy(),
                                               "offtarget_specific_parameters": off_activity.specific_parameters,
                                               "delta_transformation_parameters": delta_params
                                           })

        qed_score = ComponentParameters(component_type=enum.QED_SCORE,
                                        name="qed_score",
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters={})
        matching_substructure = ComponentParameters(component_type=enum.MATCHING_SUBSTRUCTURE,
                                                    name="matching_substructure",
                                                    weight=1.,
                                                    smiles=[AMOXAPINE],
                                                    model_path="",
                                                    specific_parameters={})
        custom_alerts = ComponentParameters(component_type=enum.CUSTOM_ALERTS,
                                            name="custom_alerts",
                                            weight=1.,
                                            smiles=[],
                                            model_path="",
                                            specific_parameters={})
        self.sf_state = CustomProduct(
            parameters=[activity, selectivity, qed_score, matching_substructure, custom_alerts])


    def test_selectivity_function_with_scikit_and_wrapped_models_1(self):
        score: FinalSummary = self.sf_state.get_final_score(smiles=[ASPIRIN])
        self.assertAlmostEqual(score.total_score[0], 0.154, 3)

    def test_selectivity_function_with_scikit_and_wrapped_models_2(self):
        score = self.sf_state.get_final_score(smiles=[CELECOXIB])
        self.assertAlmostEqual(score.total_score[0], 0.171, 3)

