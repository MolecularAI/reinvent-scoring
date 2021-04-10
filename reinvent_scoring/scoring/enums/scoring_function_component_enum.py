

class ScoringFunctionComponentNameEnum:
    PARALLEL_ROCS_SIMILARITY = "parallel_rocs_similarity"
    SELECTIVITY = "selectivity"
    PREDICTIVE_PROPERTY = "predictive_property"
    ROCS_SIMILARITY = "rocs_similarity"
    MATCHING_SUBSTRUCTURE = "matching_substructure"
    TANIMOTO_SIMILARITY = "tanimoto_similarity"
    JACCARD_DISTANCE = "jaccard_distance"
    CUSTOM_ALERTS = "custom_alerts"
    QED_SCORE = "qed_score"
    MOLECULAR_WEIGHT = "molecular_weight"
    NUM_ROTATABLE_BONDS = "num_rotatable_bonds"
    NUM_HBD_LIPINSKI = "num_hbd_lipinski"
    NUM_HBA_LIPINSKI = "num_hba_lipinski"
    NUM_RINGS = "num_rings"
    TPSA = "tpsa"
    SLOGP = "slogp"
    TOTAL_SCORE = "total_score" # there is no actual component corresponding to this type
    REACTION_FILTERS = "reaction_filters"
    #NOTE: components below are AZ specific
    AZ_LOGD74 = "az_logd74"
    HLM_CLINT = "hlm_clint"
    RH_CLINT = "rh_clint"
    HH_CLINT = "hh_clint"
    SOLUBILITY_DD = "solubilityDD"
    CACO2_INTR = "caco2_intr"
    CACO2_EFFLUX = "caco2_efflux"
    HERG = "clab_herg"
    SA_SCORE = "sa_score"
    AZDOCK = "azdock"
    AZ_LOGD74_PIP = "azlogd74"
    CACO2_INTR_PIP = "intrcaco2"
    CACO2_EFFLUX_PIP = "caco2-efflux"
    HH_CLINT_PIP = "hh-clint"
    HLM_CLINT_PIP = "hlm-clint"
    RH_CLINT_PIP = "rh-clint"
    SOLUBILITY_DD_PIP = "solubility-dd"
    HERG_PIP = "herg"
    KPUU_PIP = "kpuu-brain"
    RAT_PK_PIP = "rat-pk"
    CLAB_TOP_20 = "clab_top_20"
    RA_SCORE = "rascore"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    def __setattr__(self, key, value):
        raise ValueError("Do not assign value to a ScoringFunctionComponentNameEnum field.")
