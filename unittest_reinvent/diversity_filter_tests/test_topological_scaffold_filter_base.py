import os
import shutil
import unittest

import numpy as np

from reinvent_scoring.scoring.diversity_filters.reinvent_core.diversity_filter import DiversityFilter
from reinvent_scoring.scoring.diversity_filters.reinvent_core.diversity_filter_parameters import \
    DiversityFilterParameters
from unittest_reinvent.fixtures.paths import MAIN_TEST_PATH
from unittest_reinvent.fixtures.test_data import ASPIRIN, PROPANE
from reinvent_scoring.scoring.enums.diversity_filter_enum import DiversityFilterEnum

from reinvent_scoring.scoring.enums.scoring_function_component_enum import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.score_summary import FinalSummary, ComponentSummary
from reinvent_scoring.scoring.component_parameters import ComponentParameters


class BaseTopologicalScaffoldFilter(unittest.TestCase):

    def setUp(self):
        self.scaffold_enum = DiversityFilterEnum()
        self.sf_enum = ScoringFunctionComponentNameEnum()
        self.workfolder = MAIN_TEST_PATH
        # create a scaffold filter and fill it with a few entries
        smiles = [ASPIRIN, PROPANE]
        scores = np.array([0.7, 0.5])
        valid_idx = [0, 1]
        component_parameters = ComponentParameters(component_type=self.sf_enum.TANIMOTO_SIMILARITY,
                                                   name="tanimoto_similarity",
                                                   weight=1.,
                                                   smiles=smiles,
                                                   model_path="",
                                                   specific_parameters={})
        component_score_summary = ComponentSummary(scores, component_parameters)

        final_summary = FinalSummary(scores, smiles, valid_idx, [component_score_summary])

        sf_parameters = DiversityFilterParameters(name=self.scaffold_enum.IDENTICAL_TOPOLOGICAL_SCAFFOLD, minscore=0.5,
                                                  minsimilarity=0.4, bucket_size=1)

        self.scaffold_filter = DiversityFilter(sf_parameters)
        self.scaffold_filter.update_score(final_summary)

    def tearDown(self):
        if os.path.isdir(self.workfolder):
            shutil.rmtree(self.workfolder)
