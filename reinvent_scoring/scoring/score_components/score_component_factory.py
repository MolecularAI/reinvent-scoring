from typing import List

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import TanimotoSimilarity, \
    JaccardDistance, CustomAlerts, QedScore, MatchingSubstructure, \
    RocsSimilarity, ParallelRocsSimilarity, PredictivePropertyComponent, SelectivityComponent, \
    SASComponent, MolWeight, PSA, RotatableBonds, HBD_Lipinski, HBA_Lipinski, \
    NumRings, AZlogD74, HLMClint, SlogP, \
    RHClint, HHClint, SolubilityDD, HERG, CACO2Intrinsic, CACO2Efflux, AZdock, RatPKPiP, Top20
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.score_components.pip.pip_log_prediction_component import PiPLogPredictionComponent
from reinvent_scoring.scoring.score_components.pip.pip_prediction_component import PiPPredictionComponent


class ScoreComponentFactory:
    def __init__(self, parameters: List[ComponentParameters]):
        self._parameters = parameters
        self._current_components = self._deafult_scoring_component_registry()

    def _deafult_scoring_component_registry(self) -> dict:
        enum = ScoringFunctionComponentNameEnum()
        component_map = {
            enum.MATCHING_SUBSTRUCTURE: MatchingSubstructure,
            enum.ROCS_SIMILARITY: RocsSimilarity,
            enum.PREDICTIVE_PROPERTY: PredictivePropertyComponent,
            enum.TANIMOTO_SIMILARITY: TanimotoSimilarity,
            enum.JACCARD_DISTANCE: JaccardDistance,
            enum.CUSTOM_ALERTS: CustomAlerts,
            enum.QED_SCORE: QedScore,
            enum.MOLECULAR_WEIGHT: MolWeight,
            enum.TPSA: PSA,
            enum.NUM_ROTATABLE_BONDS: RotatableBonds,
            enum.NUM_HBD_LIPINSKI: HBD_Lipinski,
            enum.NUM_HBA_LIPINSKI: HBA_Lipinski,
            enum.NUM_RINGS: NumRings,
            enum.SLOGP: SlogP,
            enum.PARALLEL_ROCS_SIMILARITY: ParallelRocsSimilarity,
            enum.AZ_LOGD74: AZlogD74,
            enum.HLM_CLINT: HLMClint,
            enum.RH_CLINT: RHClint,
            enum.HH_CLINT: HHClint,
            enum.SOLUBILITY_DD: SolubilityDD,
            enum.HERG: HERG,
            enum.CACO2_INTR: CACO2Intrinsic,
            enum.CACO2_EFFLUX: CACO2Efflux,
            enum.SELECTIVITY: SelectivityComponent,
            enum.SA_SCORE: SASComponent,
            enum.AZDOCK: AZdock,
            enum.AZ_LOGD74_PIP: PiPPredictionComponent,
            enum.CACO2_INTR_PIP: PiPPredictionComponent,
            enum.CACO2_EFFLUX_PIP: PiPLogPredictionComponent,
            enum.HH_CLINT_PIP: PiPLogPredictionComponent,
            enum.HLM_CLINT_PIP: PiPLogPredictionComponent,
            enum.RH_CLINT_PIP: PiPLogPredictionComponent,
            enum.SOLUBILITY_DD_PIP: PiPLogPredictionComponent,
            enum.HERG_PIP: PiPPredictionComponent,
            enum.RAT_PK_PIP: RatPKPiP,
            enum.CLAB_TOP_20: Top20,
            enum.RA_SCORE: PiPPredictionComponent,
            enum.KPUU_PIP: PiPLogPredictionComponent
        }
        return component_map

    def create_score_components(self) -> [BaseScoreComponent]:
        return [self._current_components.get(p.component_type)(p) for p in self._parameters]
