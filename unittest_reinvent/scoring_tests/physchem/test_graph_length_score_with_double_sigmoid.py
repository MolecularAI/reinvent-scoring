import numpy as np
import numpy.testing as npt

from unittest_reinvent.fixtures.test_data import CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE
from unittest_reinvent.scoring_tests.physchem.base_setup import BaseSetup


class TestGraphLengthScoreWithDoubleSigmoid(BaseSetup):

    def setUp(self):
        super().setup_attrs()
        specific_parameters = {
            self.csp_enum.TRANSFORMATION: True,
            self.csp_enum.LOW: 10,
            self.csp_enum.HIGH: 20,
            self.csp_enum.COEF_DIV: 25,
            self.csp_enum.COEF_SI: 20,
            self.csp_enum.COEF_SE: 20,
            self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.DOUBLE_SIGMOID
        }
        super().init(self.sf_enum.GRAPH_LENGTH, specific_parameters)
        super().setUp()

    def test_graph_length_1(self):
        smiles = [CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE]
        values = np.array([0.975, 0.5, 0.137, 0, 0.5])
        score = self.sf_state.get_final_score(smiles=smiles)
        npt.assert_array_almost_equal(score.total_score, values, 2)
