import numpy as np
import numpy.testing as npt

from unittest_reinvent.fixtures.test_data import CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE
from unittest_reinvent.scoring_tests.physchem.base_setup import BaseSetup


class TestSlogpScoreWithDoubleSigmoid(BaseSetup):

    def setUp(self):
        super().setup_attrs()
        specific_parameters = {
            self.csp_enum.TRANSFORMATION: True,
            self.csp_enum.LOW: 1,
            self.csp_enum.HIGH: 4,
            self.csp_enum.COEF_DIV: 3,
            self.csp_enum.COEF_SI: 10,
            self.csp_enum.COEF_SE: 10,
            self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.DOUBLE_SIGMOID
        }
        super().init(self.sf_enum.SLOGP, specific_parameters)
        super().setUp()

    def test_slogp_1(self):
        smiles = [CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE]
        values = np.array([0.977, 0.142, 0.988, 0, 1])
        score = self.sf_state.get_final_score(smiles=smiles)
        npt.assert_array_almost_equal(score.total_score, values, 2)
