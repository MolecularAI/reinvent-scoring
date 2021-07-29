import numpy as np
import numpy.testing as npt

from reinvent_scoring import TransformationTypeEnum, ScoringFunctionComponentNameEnum
from unittest_reinvent.fixtures.test_data import CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE
from unittest_reinvent.scoring_tests.physchem.base_setup import BaseSetup


class TestHbdScoreWithDoubleSigmoid(BaseSetup):

    def setUp(self):
        super().setup_attrs()
        specific_parameters = {
            self.csp_enum.TRANSFORMATION: True,
            self.csp_enum.LOW: 0,
            self.csp_enum.HIGH: 1,
            self.csp_enum.TRANSFORMATION_TYPE: self.tt_enum.STEP
        }
        super().init(self.sf_enum.NUM_HBD_LIPINSKI, specific_parameters)
        super().setUp()

    def test_hbd_1(self):
        smiles = [CELECOXIB, METAMIZOLE, AMOXAPINE, METHOXYHYDRAZINE, COCAINE]
        values = np.array([1., 1., 1., 0., 1.])
        score = self.sf_state.get_final_score(smiles=smiles)
        npt.assert_array_equal(score.total_score, values)
