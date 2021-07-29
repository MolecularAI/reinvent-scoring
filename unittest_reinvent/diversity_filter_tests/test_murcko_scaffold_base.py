import unittest

from reinvent_scoring.scoring.diversity_filters.reinvent_core.diversity_filter import DiversityFilter
from reinvent_scoring.scoring.diversity_filters.reinvent_core.diversity_filter_parameters import \
    DiversityFilterParameters
from unittest_reinvent.fixtures.test_data import PROPANE, PENTANE, ASPIRIN


from unittest_reinvent.diversity_filter_tests.fixtures import tanimoto_scaffold_filter_arrangement
from reinvent_scoring.scoring.enums.diversity_filter_enum import DiversityFilterEnum

from reinvent_scoring.scoring.enums.scoring_function_component_enum import ScoringFunctionComponentNameEnum


class BaseMurckoScaffoldFilter(unittest.TestCase):
    """Base class that other classes in test_murcko_scaffold files inherit from."""

    def setUp(self):
        self.scaffold_enum = DiversityFilterEnum()
        self.sf_enum = ScoringFunctionComponentNameEnum()

        final_summary = tanimoto_scaffold_filter_arrangement([ASPIRIN, PROPANE, PENTANE], [0.7, 0.5, 0.], [0, 1, 2])

        sf_parameters = DiversityFilterParameters(name=self.scaffold_enum.IDENTICAL_MURCKO_SCAFFOLD, minscore=0.5,
                                                  minsimilarity=0.4, bucket_size=1)

        self.scaffold_filter = DiversityFilter(sf_parameters)
        self.scaffold_filter.update_score(final_summary)
