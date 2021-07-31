from reinvent_scoring import TransformationFactory
from unittest_reinvent.scoring_tests.scoring_transformations.base_test_score_transformations import BaseTestScoreTransformation


class BaseTestSigmoidScoreTransformation(BaseTestScoreTransformation):

    def setUp(self):
        self.v_list = [12.016, 0.0015, 7.9, 123.264, 77.80, 4.0, 111.12]
        self.factory = TransformationFactory()
        base_specific_parameters = {self.csp_enum.TRANSFORMATION: True,
                                    self.csp_enum.LOW: 4,
                                    self.csp_enum.HIGH: 9,
                                    self.csp_enum.K: 0.25}

        self.specific_parameters = {**base_specific_parameters, **self.added_specific_parameters}

        transform_function = self.factory.get_transformation_function(self.specific_parameters)
        self.transformed_scores = transform_function(predictions=self.v_list[:],
                                                     parameters=self.specific_parameters)

