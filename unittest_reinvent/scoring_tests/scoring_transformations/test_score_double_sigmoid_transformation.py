from unittest_reinvent.scoring_tests.scoring_transformations.base_test_score_sigmoid_transformations import \
    BaseTestSigmoidScoreTransformation
from unittest_reinvent.scoring_tests.scoring_transformations.fixtures import round_list


class TestDoubleSigmoidTransformations(BaseTestSigmoidScoreTransformation):

    def setUp(self):
        super().setup_attrs()
        specific_parameters = {self.csp_enum.COEF_DIV: 100,
                               self.csp_enum.COEF_SI: 150,
                               self.csp_enum.COEF_SE: 150,
                               self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.DOUBLE_SIGMOID}
        super().init(specific_parameters)
        super().setUp()

    def test_double_sigmoid_transformation(self):
        self.assertListEqual(
            [0, 0, 0.978, 0.0, 0.0, 0.5, 0.0],
            round_list(self.transformed_scores.tolist(), 3))

    def test_double_sigmoid_transformation_updated_parameters(self):
        self.update_parameters({self.csp_enum.COEF_DIV: 200, self.csp_enum.COEF_SI: 100, self.csp_enum.COEF_SE: 100})
        transform_function = self.factory.get_transformation_function(self.specific_parameters)
        self.transformed_scores = transform_function(predictions=self.v_list[:],
                                                     parameters=self.specific_parameters)
        self.assertListEqual(
            [0.03, 0.01, 0.769, 0.0, 0.0, 0.497, 0.0],
            round_list(self.transformed_scores.tolist(), 3))
