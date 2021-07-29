from unittest_reinvent.diversity_filter_tests.test_tanimoto_similarity_base import BaseTanimotoSimilarity
from unittest_reinvent.fixtures.test_data import CELECOXIB
from unittest_reinvent.diversity_filter_tests.fixtures import tanimoto_scaffold_filter_arrangement


class TestScaffoldSimilarityFilterValidAddition(BaseTanimotoSimilarity):

    def setUp(self):
        super().setUp()

        self.final_summary = tanimoto_scaffold_filter_arrangement([CELECOXIB], [1.0], [0])

    def test_valid_addition(self):
        self.scaffold_filter.update_score(self.final_summary)
        self.assertEqual(3, self.scaffold_filter._diversity_filter_memory.number_of_scaffolds())
