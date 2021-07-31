from unittest_reinvent.scoring_tests.scoring_transformations.base_test_score_sigmoid_transformations import \
    BaseTestSigmoidScoreTransformation
from unittest_reinvent.scoring_tests.scoring_transformations.fixtures import round_list


class TestReverseSigmoidTransformations(BaseTestSigmoidScoreTransformation):

    def setUp(self):
        super().setup_attrs()
        specific_parameters = {self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.REVERSE_SIGMOID}
        super().init(specific_parameters)
        super().setUp()

    def test_reverse_sigmoid_transformation(self):
        self.assertListEqual(
            [0.002, 0.999, 0.166, 0.0, 0, 0.947, 0.0],
            round_list(self.transformed_scores.tolist(), 3))

    def test_reverse_sigmoid_transformation_updated_parameters(self):
        self.update_parameters({self.csp_enum.HIGH: 45, self.csp_enum.LOW: 9})
        self.assertListEqual(
            [0.917, 0.987, 0.955, 0, 0, 0.975, 0],
            round_list(self.transformed_scores.tolist(), 3))
