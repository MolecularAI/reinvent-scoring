import unittest

from reinvent_scoring import TransformationTypeEnum, TransformationFactory
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum


class BaseTestScoreTransformation(unittest.TestCase):

    def setup_attrs(self):
        self.tt_enum = TransformationTypeEnum()
        self.csp_enum = ComponentSpecificParametersEnum()

    def init(self, added_specific_parameters):
        self.added_specific_parameters = added_specific_parameters

    def setUp(self):
        self.v_list = [12.086, 0.0015, 7.9, 123.264, 77.80, 4.0, 111.12]
        self.factory = TransformationFactory()
        base_specific_parameters = {self.csp_enum.TRANSFORMATION: True,
                                    self.csp_enum.LOW: 4}

        self.specific_parameters = {**base_specific_parameters, **self.added_specific_parameters}

        transform_function = self.factory.get_transformation_function(self.specific_parameters)
        self.transformed_scores = transform_function(predictions=self.v_list[:],
                                                     parameters=self.specific_parameters)

    def update_parameters(self, parameters):
        for parameter in parameters:
            self.specific_parameters[parameter] = parameters[parameter]

        transform_function = self.factory.get_transformation_function(self.specific_parameters)
        self.transformed_scores = transform_function(predictions=self.v_list[:],
                                                     parameters=self.specific_parameters)
