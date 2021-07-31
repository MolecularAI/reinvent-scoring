import unittest
from typing import List

from reinvent_scoring import ComponentSpecificParametersEnum, TransformationTypeEnum, ScoringFunctionComponentNameEnum, \
    ComponentParameters
from unittest_reinvent.scoring_tests.fixtures import create_activity_component_regression, \
    create_offtarget_activity_component_regression


class BaseTestSelectivityFunctionDoubleSigmoid(unittest.TestCase):

    def init(self, selectivity_smiles: List[str] = [], qed_smiles: List[str] = [],
             matching_substructure_smiles: List[str] = [], custom_alert_smiles: List[str] = [],):
        self.selectivity_smiles = selectivity_smiles
        self.qed_smiles = qed_smiles
        self.matching_substructure_smiles = matching_substructure_smiles
        self.custom_alert_smiles = custom_alert_smiles

    def setUp(self):
        self.csp_enum = ComponentSpecificParametersEnum()
        self.transf_type = TransformationTypeEnum()
        self.sf_enum = ScoringFunctionComponentNameEnum()
        self.activity = create_activity_component_regression()
        self.activity.specific_parameters[self.csp_enum.TRANSFORMATION_TYPE] = self.transf_type.DOUBLE_SIGMOID
        self.activity.specific_parameters[self.csp_enum.COEF_DIV] = 100.
        self.activity.specific_parameters[self.csp_enum.COEF_SI] = 150.
        self.activity.specific_parameters[self.csp_enum.COEF_SE] = 150.
        self.off_activity = create_offtarget_activity_component_regression()

        self.delta_params = {
            "high": 3.0,
            "k": 0.25,
            "low": 0.0,
            "transformation": True,
            "transformation_type": "sigmoid"
        }

        self.selectivity = ComponentParameters(component_type=self.sf_enum.SELECTIVITY,
                                               name="desirability",
                                               weight=1.,
                                               smiles=self.selectivity_smiles,
                                               model_path="",
                                               specific_parameters={
                                                   "activity_model_path": self.activity.model_path,
                                                   "offtarget_model_path": self.off_activity.model_path,
                                                   "activity_specific_parameters": self.activity.specific_parameters.copy(),
                                                   "offtarget_specific_parameters": self.off_activity.specific_parameters,
                                                   "delta_transformation_parameters": self.delta_params
                                               })

        self.qed_score = ComponentParameters(component_type=self.sf_enum.QED_SCORE,
                                             name="qed_score",
                                             weight=1.,
                                             smiles=self.qed_smiles,
                                             model_path="",
                                             specific_parameters={})
        self.matching_substructure = ComponentParameters(component_type=self.sf_enum.MATCHING_SUBSTRUCTURE,
                                                         name="matching_substructure",
                                                         weight=1.,
                                                         smiles=self.matching_substructure_smiles,
                                                         model_path="",
                                                         specific_parameters={})
        self.custom_alerts = ComponentParameters(component_type=self.sf_enum.CUSTOM_ALERTS,
                                                 name="custom_alerts",
                                                 weight=1.,
                                                 smiles=self.custom_alert_smiles,
                                                 model_path="",
                                                 specific_parameters={})
