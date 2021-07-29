from unittest_reinvent.diversity_filter_tests.test_tanimoto_similarity_base import BaseTanimotoSimilarity
from unittest_reinvent.fixtures.test_data import INVALID
from unittest_reinvent.diversity_filter_tests.fixtures import tanimoto_scaffold_filter_arrangement


class TestNoFilterInvalidSmileProtection(BaseTanimotoSimilarity):

    def setUp(self):
        super().setUp()

        self.final_summary = tanimoto_scaffold_filter_arrangement([INVALID], [1.0], [0])

    def test_superfluous_addition(self):
        try:
            self.scaffold_filter.update_score(self.final_summary)
        except Exception as e:
            self.assertEqual(type(e).__name__, "ArgumentError")
        else:
            self.fail("""Expected exception of type "ArgumentError" because of invalid smile.""")