class ROCSSpecificParametersEnum():
    _ROCS_INPUT = "rocs_input"
    _INPUT_TYPE = "input_type"
    _SHAPE_WEIGHT = "shape_weight"
    _COLOR_WEIGHT = "color_weight"
    _SIM_MEASURE = "similarity_measure"
    _MAX_CPUS = "max_num_cpus"
    _CUSTOM_CFF = "custom_cff"
    _SAVE_ROCS_OVERLAYS = "save_rocs_overlays"
    _ROCS_OVERLAYS_DIR = "rocs_overlays_dir"
    _ROCS_OVERLAYS_PREFIX = "rocs_overlays_prefix"
    _ENUM_STEREO = "enumerate_stereo"
    _MAX_STEREO = "max_stereocenters"
    _NEGATIVE_VOLUME = "negative_volume"
    _PROTEIN_NEG_VOL_FILE = "protein_neg_vol_file"
    _LIGAND_NEG_VOL_FILE = "ligand_neg_vol_file"
    _MAX_CONFS = "max_confs"
    _EWINDOW = "ewindow"

    # Need to stick to this pattern of getter/setter as the alternative more concise one cannot
    # be pickled for multiprocessing

    @property
    def ROCS_INPUT(self):
        return self._ROCS_INPUT

    @ROCS_INPUT.setter
    def ROCS_INPUT(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def INPUT_TYPE(self):
        return self._INPUT_TYPE

    @INPUT_TYPE.setter
    def INPUT_TYPE(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def SHAPE_WEIGHT(self):
        return self._SHAPE_WEIGHT

    @SHAPE_WEIGHT.setter
    def SHAPE_WEIGHT(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def COLOR_WEIGHT(self):
        return self._COLOR_WEIGHT

    @COLOR_WEIGHT.setter
    def COLOR_WEIGHT(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def SIM_MEASURE(self):
        return self._SIM_MEASURE

    @SIM_MEASURE.setter
    def SIM_MEASURE(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def MAX_CPUS(self):
        return self._MAX_CPUS

    @MAX_CPUS.setter
    def MAX_CPUS(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def CUSTOM_CFF(self):
        return self._CUSTOM_CFF

    @CUSTOM_CFF.setter
    def CUSTOM_CFF(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def SAVE_ROCS_OVERLAYS(self):
        return self._SAVE_ROCS_OVERLAYS

    @SAVE_ROCS_OVERLAYS.setter
    def SAVE_ROCS_OVERLAYS(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def ROCS_OVERLAYS_DIR(self):
        return self._ROCS_OVERLAYS_DIR

    @ROCS_OVERLAYS_DIR.setter
    def ROCS_OVERLAYS_DIR(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def ROCS_OVERLAYS_PREFIX(self):
        return self._ROCS_OVERLAYS_PREFIX

    @ROCS_OVERLAYS_PREFIX.setter
    def ROCS_OVERLAYS_PREFIX(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def ENUM_STEREO(self):
        return self._ENUM_STEREO

    @ENUM_STEREO.setter
    def ENUM_STEREO(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def MAX_STEREO(self):
        return self._MAX_STEREO

    @MAX_STEREO.setter
    def MAX_STEREO(self, value):
        raise ValueError("Do not assign value to a ROCSSpecificParametersEnum field")

    @property
    def PROTEIN_NEG_VOL_FILE(self):
        return self._PROTEIN_NEG_VOL_FILE

    @PROTEIN_NEG_VOL_FILE.setter
    def PROTEIN_NEG_VOL_FILE(self, value):
        raise ValueError("Do not assign value to a ROCSInputFileTypesEnum field")

    @property
    def NEGATIVE_VOLUME(self):
        return self._NEGATIVE_VOLUME

    @NEGATIVE_VOLUME.setter
    def NEGATIVE_VOLUME(self, value):
        raise ValueError("Do not assign value to a ROCSInputFileTypesEnum field")

    @property
    def LIGAND_NEG_VOL_FILE(self):
        return self._LIGAND_NEG_VOL_FILE

    @LIGAND_NEG_VOL_FILE.setter
    def LIGAND_NEG_VOL_FILE(self, value):
        raise ValueError("Do not assign value to a ROCSInputFileTypesEnum field")

    @property
    def MAX_CONFS(self):
        return self._MAX_CONFS

    @MAX_CONFS.setter
    def MAX_CONFS(self, value):
        raise ValueError("Do not assign value to a ROCSInputFileTypesEnum field")

    @property
    def EWINDOW(self):
        return self._EWINDOW

    @EWINDOW.setter
    def EWINDOW(self, value):
        raise ValueError("Do not assign value to a ROCSInputFileTypesEnum field")