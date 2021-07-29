import unittest

from reinvent_scoring import ScoringFunctionComponentNameEnum, ComponentSpecificParametersEnum, TransformationTypeEnum, \
    ComponentParameters, SelectivityComponent


class BaseTestSelectivityComponent(unittest.TestCase):
    def setUp(self, broken: bool = False):
        csp_enum = ComponentSpecificParametersEnum()
        transf_type = TransformationTypeEnum()
        enum = ScoringFunctionComponentNameEnum()

        delta_params = {
            "high": 3.0,
            "k": 0.25,
            "low": 0.0,
            "transformation": True,
            "transformation_type": "sigmoid"
        }
        activity = self.activity
        activity.specific_parameters[csp_enum.TRANSFORMATION_TYPE] = transf_type.DOUBLE_SIGMOID
        activity.specific_parameters[csp_enum.COEF_DIV] = 100.
        activity.specific_parameters[csp_enum.COEF_SI] = 150.
        activity.specific_parameters[csp_enum.COEF_SE] = 150.

        off_activity = self.off_activity

        if broken:
            activity.specific_parameters.pop(csp_enum.TRANSFORMATION_TYPE, None)
            off_activity.specific_parameters.pop(csp_enum.TRANSFORMATION_TYPE, None)

        selectivity = ComponentParameters(component_type=enum.SELECTIVITY,
                                          name="desirability",
                                          weight=1.,
                                          smiles=[],
                                          model_path="",
                                          specific_parameters={
                                               "activity_model_path": activity.model_path,
                                               "offtarget_model_path": off_activity.model_path,
                                               "activity_specific_parameters": activity.specific_parameters.copy(),
                                               "offtarget_specific_parameters": off_activity.specific_parameters.copy(),
                                               "delta_transformation_parameters": delta_params
                                           })

        self.component = SelectivityComponent(parameters=selectivity)

