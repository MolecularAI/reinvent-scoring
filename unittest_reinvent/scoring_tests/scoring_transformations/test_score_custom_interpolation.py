import unittest

from reinvent_scoring import TransformationTypeEnum, ComponentSpecificParametersEnum, TransformationFactory


class TestScoreCustomInterpolation(unittest.TestCase):

    def setUp(self):
        self.tt_enum = TransformationTypeEnum()
        self.csp_enum = ComponentSpecificParametersEnum()
        specific_parameters = {self.csp_enum.TRUNCATE_RIGHT: True,
                               self.csp_enum.TRUNCATE_LEFT: True,
                               self.csp_enum.INTERPOLATION_MAP: [{"origin": 0.0, "destination": 0.0},
                                                                 {"origin": 1.0, "destination": 1.0}],
                               self.csp_enum.TRANSFORMATION: True,
                               self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.CUSTOM_INTERPOLATION}
        self.factory = TransformationFactory()
        transform_function = self.factory.get_transformation_function(specific_parameters)
        self.transformed_scores = transform_function(predictions=[-1, 0., 0.3, 0.7, 15],
                                                     parameters=specific_parameters)

    def test_custom_interpolation(self):
        self.assertListEqual([0, 0, 0.3, 0.7, 1], self.transformed_scores.tolist())
