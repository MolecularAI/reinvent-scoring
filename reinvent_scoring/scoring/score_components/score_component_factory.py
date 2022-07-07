from typing import List

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_components import TanimotoSimilarity, \
    JaccardDistance, CustomAlerts, QedScore, MatchingSubstructure, \
    RocsSimilarity, ParallelRocsSimilarity, PredictivePropertyComponent, SelectivityComponent, \
    SASComponent, MolWeight, PSA, RotatableBonds, HBD_Lipinski, HBA_Lipinski, \
    NumRings, SlogP, AZdock, RatPKPiP, PiPLogPredictionComponent, PiPPredictionComponent, \
    QptunaPiPModelComponent, StringPiPPredictionComponent, GraphLength, NumberOfStereoCenters, \
    LinkerLengthRatio, LinkerGraphLength, LinkerEffectiveLength, LinkerNumRings, LinkerNumAliphaticRings, \
    LinkerNumAromaticRings, LinkerNumSPAtoms, LinkerNumSP2Atoms, LinkerNumSP3Atoms, LinkerNumHBA, \
    LinkerNumHBD, LinkerMolWeight, LinkerRatioRotatableBonds, DockStream, NumAromaticRings, NumAliphaticRings
from reinvent_scoring.scoring.score_components.console_invoked import Icolos
from reinvent_scoring.scoring.score_components.dlf.dlf_api_score_component import DlfApiComponent


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
            enum.DLF: DlfApiComponent,
            enum.JACCARD_DISTANCE: JaccardDistance,
            enum.CUSTOM_ALERTS: CustomAlerts,
            enum.QED_SCORE: QedScore,
            enum.MOLECULAR_WEIGHT: MolWeight,
            enum.TPSA: PSA,
            enum.NUM_ROTATABLE_BONDS: RotatableBonds,
            enum.GRAPH_LENGTH: GraphLength,
            enum.NUM_HBD_LIPINSKI: HBD_Lipinski,
            enum.NUM_HBA_LIPINSKI: HBA_Lipinski,
            enum.NUM_RINGS: NumRings,
            enum.NUM_AROMATIC_RINGS: NumAromaticRings,
            enum.NUM_ALIPHATIC_RINGS: NumAliphaticRings,
            enum.SLOGP: SlogP,
            enum.NUMBER_OF_STEREO_CENTERS: NumberOfStereoCenters,
            enum.PARALLEL_ROCS_SIMILARITY: ParallelRocsSimilarity,
            enum.SELECTIVITY: SelectivityComponent,
            enum.SA_SCORE: SASComponent,
            enum.AZDOCK: AZdock,
            enum.AZ_LOGD74_PIP: PiPPredictionComponent,
            enum.CACO2_INTR_PIP: PiPLogPredictionComponent,
            enum.CACO2_EFFLUX_PIP: PiPPredictionComponent,
            enum.HH_CLINT_PIP: PiPLogPredictionComponent,
            enum.HLM_CLINT_PIP: PiPLogPredictionComponent,
            enum.RH_CLINT_PIP: PiPLogPredictionComponent,
            enum.SOLUBILITY_DD_PIP: PiPLogPredictionComponent,
            enum.HERG_PIP: PiPPredictionComponent,
            enum.RAT_PK_PIP: RatPKPiP,
            enum.RA_SCORE: PiPPredictionComponent,
            enum.KPUU_PIP: PiPLogPredictionComponent,
            enum.QPTUNA_PIP_MODEL: QptunaPiPModelComponent,
            enum.THP1_CYTOTOXICITY: StringPiPPredictionComponent,
            enum.LINKER_GRAPH_LENGTH: LinkerGraphLength,
            enum.LINKER_EFFECTIVE_LENGTH: LinkerEffectiveLength,
            enum.LINKER_LENGTH_RATIO: LinkerLengthRatio,
            enum.LINKER_NUM_RINGS: LinkerNumRings,
            enum.LINKER_NUM_ALIPHATIC_RINGS: LinkerNumAliphaticRings,
            enum.LINKER_NUM_AROMATIC_RINGS: LinkerNumAromaticRings,
            enum.LINKER_NUM_SP_ATOMS: LinkerNumSPAtoms,
            enum.LINKER_NUM_SP2_ATOMS: LinkerNumSP2Atoms,
            enum.LINKER_NUM_SP3_ATOMS: LinkerNumSP3Atoms,
            enum.LINKER_NUM_HBA: LinkerNumHBA,
            enum.LINKER_NUM_HBD: LinkerNumHBD,
            enum.LINKER_MOL_WEIGHT: LinkerMolWeight,
            enum.LINKER_RATIO_ROTATABLE_BONDS: LinkerRatioRotatableBonds,
            enum.DOCKSTREAM: DockStream,
            enum.ICOLOS: Icolos,
            # enum.AIZYNTH: BuildingBlockAvailabilityComponent
        }
        return component_map

    def create_score_components(self) -> [BaseScoreComponent]:
        def create_component(component_params):
            if component_params.component_type in self._current_components:
                component = self._current_components[component_params.component_type]
                component_instance = component(component_params)
            else:
                raise KeyError(f'Component: {component_params.component_type} is not implemented.'
                               f' Consider checking your input.')
            return component_instance

        components = [create_component(component) for component in self._parameters]
        return components
