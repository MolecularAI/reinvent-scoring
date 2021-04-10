class ContainerType:
    _SCIKIT_CONTAINER = "scikit_container"
    _OPTUNA_CONTAINER = "optuna_container"

    @property
    def SCIKIT_CONTAINER(self):
        return self._SCIKIT_CONTAINER

    @SCIKIT_CONTAINER.setter
    def SCIKIT_CONTAINER(self, value):
        raise ValueError("Do not assign value to a ContainerType field")

    @property
    def OPTUNA_CONTAINER(self):
        return self._OPTUNA_CONTAINER

    @OPTUNA_CONTAINER.setter
    def OPTUNA_CONTAINER(self, value):
        raise ValueError("Do not assign value to a ContainerType field")
