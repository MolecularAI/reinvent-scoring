import numpy as np
import numpy.testing as npt

from unittest_reinvent.fixtures.test_data import CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE
from unittest_reinvent.scoring_tests.physchem.base_setup import BaseSetup


class TestMwScoreWithDoubleSigmoid(BaseSetup):

    def setUp(self):
        super().setup_attrs()
        specific_parameters = {
            self.csp_enum.TRANSFORMATION: True,
            self.csp_enum.LOW: 200,
            self.csp_enum.HIGH: 400,
            self.csp_enum.COEF_DIV: 500,
            self.csp_enum.COEF_SI: 20,
            self.csp_enum.COEF_SE: 20,
            self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.DOUBLE_SIGMOID
        }
        super().init(self.sf_enum.MOLECULAR_WEIGHT, specific_parameters)
        super().setUp()

    def test_mw_1(self):
        smiles = [CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE]
        values = np.array([8.47e-01, 1.00e+00, 1.00e+00, 3.04e-06, 1.00e+00])
        score = self.sf_state.get_final_score(smiles=smiles)
        npt.assert_array_almost_equal(score.total_score, values, 2)
