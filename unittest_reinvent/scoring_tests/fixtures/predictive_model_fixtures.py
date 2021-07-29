from reinvent_scoring.scoring.component_parameters import ComponentParameters
from unittest_reinvent.fixtures.paths import ACTIVITY_REGRESSION, ACTIVITY_CLASSIFICATION, MAIN_TEST_PATH
from unittest_reinvent.fixtures.paths import AZDOCK_DOCKER_SCRIPT_PATH, AZDOCK_ENV_PATH, AZDOCK_DEBUG
from unittest_reinvent.fixtures.paths import DOCKSTREAM_DOCKER_SCRIPT_PATH, DOCKSTREAM_ENV_PATH, DOCKSTREAM_DEBUG
from unittest_reinvent.fixtures.paths import AZGARD_DEBUG, AZGARD_EXECUTOR_SCRIPT_PATH, AZGARD_ENV_PATH
from reinvent_scoring.scoring.enums import ComponentSpecificParametersEnum
from reinvent_scoring.scoring.enums import DescriptorTypesEnum
from reinvent_scoring.scoring.enums import ScoringFunctionComponentNameEnum
from reinvent_scoring.scoring.enums import TransformationTypeEnum


def create_activity_component_regression():
    sf_enum = ScoringFunctionComponentNameEnum()
    specific_parameters = _specific_parameters_regression_predictive_model()
    parameters = ComponentParameters(component_type=sf_enum.PREDICTIVE_PROPERTY,
                                        name="activity",
                                        weight=1.,
                                        smiles=[],
                                        model_path=ACTIVITY_REGRESSION,
                                        specific_parameters=specific_parameters)
    return parameters


def create_AZdock_component_parameters():
    sf_enum = ScoringFunctionComponentNameEnum()
    parameters = ComponentParameters(component_type=sf_enum.AZDOCK,
                                     name=sf_enum.AZDOCK,
                                     weight=1.,
                                     smiles=[],
                                     model_path="",
                                     specific_parameters=_specific_parameters_AZdock())
    return parameters


def create_DockStream_component_parameters():
    sf_enum = ScoringFunctionComponentNameEnum()
    parameters = ComponentParameters(component_type=sf_enum.DOCKSTREAM,
                                     name=sf_enum.DOCKSTREAM,
                                     weight=1.,
                                     smiles=[],
                                     model_path="",
                                     specific_parameters=_specific_parameters_DockStream())
    return parameters


def create_AZgard_component_parameters():
    sf_enum = ScoringFunctionComponentNameEnum()
    parameters = ComponentParameters(component_type=sf_enum.AZGARD,
                                     name=sf_enum.AZGARD,
                                     weight=1.,
                                     smiles=[],
                                     model_path="",
                                     specific_parameters=_specific_parameters_AZgard())
    return parameters


def create_offtarget_activity_component_regression():
    csp_enum = ComponentSpecificParametersEnum()
    sf_enum = ScoringFunctionComponentNameEnum()
    specific_parameters = _specific_parameters_regression_predictive_model()
    specific_parameters[csp_enum.HIGH] = 3
    specific_parameters[csp_enum.LOW] = 0
    specific_parameters[csp_enum.TRANSFORMATION] = False
    parameters = ComponentParameters(component_type=sf_enum.PREDICTIVE_PROPERTY,
                                       name="offtarget_activity",
                                       weight=1.,
                                       smiles=[],
                                       model_path=ACTIVITY_REGRESSION,
                                       specific_parameters=specific_parameters)
    return parameters


def create_predictive_property_component_regression():
    sf_enum = ScoringFunctionComponentNameEnum()
    specific_parameters = _specific_parameters_regression_predictive_model()
    parameters = ComponentParameters(component_type=sf_enum.PREDICTIVE_PROPERTY,
                                       name="predictive_property",
                                       weight=1.,
                                       smiles=[],
                                       model_path=ACTIVITY_REGRESSION,
                                       specific_parameters=specific_parameters)
    return parameters


def create_activity_component_classification():
    sf_enum = ScoringFunctionComponentNameEnum()
    specific_parameters = _specific_parameters_classifiaction_predictive_model()
    parameters = ComponentParameters(component_type=sf_enum.PREDICTIVE_PROPERTY,
                                        name="activity_classification",
                                        weight=1.,
                                        smiles=[],
                                        model_path=ACTIVITY_CLASSIFICATION,
                                        specific_parameters=specific_parameters)
    return parameters


def create_offtarget_activity_component_classification():
    csp_enum = ComponentSpecificParametersEnum()
    sf_enum = ScoringFunctionComponentNameEnum()
    specific_parameters = _specific_parameters_classifiaction_predictive_model()
    specific_parameters[csp_enum.HIGH] = 3
    specific_parameters[csp_enum.LOW] = 0
    specific_parameters[csp_enum.TRANSFORMATION] = False
    parameters = ComponentParameters(component_type=sf_enum.PREDICTIVE_PROPERTY,
                                        name="predictive_property_classification",
                                        weight=1.,
                                        smiles=[],
                                        model_path=ACTIVITY_CLASSIFICATION,
                                        specific_parameters=specific_parameters)
    return parameters


def create_predictive_property_component_classification():
    sf_enum = ScoringFunctionComponentNameEnum()
    specific_parameters = _specific_parameters_classifiaction_predictive_model()
    parameters = ComponentParameters(component_type=sf_enum.PREDICTIVE_PROPERTY,
                                        name="predictive_property_classification",
                                        weight=1.,
                                        smiles=[],
                                        model_path=ACTIVITY_CLASSIFICATION,
                                        specific_parameters=specific_parameters)
    return parameters


def create_c_lab_component(somponent_type):
    csp_enum = ComponentSpecificParametersEnum()
    transf_type = TransformationTypeEnum()
    specific_parameters = {csp_enum.CLAB_INPUT_FILE: f"{MAIN_TEST_PATH}/clab_input.json",
                           csp_enum.HIGH: 9,
                           csp_enum.LOW: 4,
                           csp_enum.K: 0.25,
                           csp_enum.TRANSFORMATION: True,
                           csp_enum.TRANSFORMATION_TYPE: transf_type.SIGMOID}
    parameters = ComponentParameters(component_type=somponent_type,
                                        name=somponent_type,
                                        weight=1.,
                                        smiles=[],
                                        model_path="",
                                        specific_parameters=specific_parameters)
    return parameters


def _specific_parameters_AZdock():
    csp_enum = ComponentSpecificParametersEnum()
    transf_type = TransformationTypeEnum()
    specific_parameters = {csp_enum.AZDOCK_CONFPATH: None,
                           csp_enum.AZDOCK_DOCKERSCRIPTPATH: AZDOCK_DOCKER_SCRIPT_PATH,
                           csp_enum.AZDOCK_ENVPATH: AZDOCK_ENV_PATH,
                           csp_enum.AZDOCK_DEBUG: AZDOCK_DEBUG,
                           csp_enum.TRANSFORMATION: True,
                           csp_enum.HIGH: -5,
                           csp_enum.LOW: -45,
                           csp_enum.K: 0.25,
                           csp_enum.TRANSFORMATION_TYPE: transf_type.REVERSE_SIGMOID}
    return specific_parameters


def _specific_parameters_DockStream():
    csp_enum = ComponentSpecificParametersEnum()
    transf_type = TransformationTypeEnum()
    specific_parameters = {csp_enum.DOCKSTREAM_CONFPATH: None,
                           csp_enum.DOCKSTREAM_DOCKERSCRIPTPATH: DOCKSTREAM_DOCKER_SCRIPT_PATH,
                           csp_enum.DOCKSTREAM_ENVPATH: DOCKSTREAM_ENV_PATH,
                           csp_enum.DOCKSTREAM_DEBUG: DOCKSTREAM_DEBUG,
                           csp_enum.TRANSFORMATION: True,
                           csp_enum.HIGH: -5,
                           csp_enum.LOW: -45,
                           csp_enum.K: 0.25,
                           csp_enum.TRANSFORMATION_TYPE: transf_type.REVERSE_SIGMOID}
    return specific_parameters


def _specific_parameters_AZgard():
    csp_enum = ComponentSpecificParametersEnum()
    specific_parameters = {csp_enum.AZGARD_CONFPATH: None,
                           csp_enum.AZGARD_EXECUTORSCRIPTPATH: AZGARD_EXECUTOR_SCRIPT_PATH,
                           csp_enum.AZGARD_ENVPATH: AZGARD_ENV_PATH,
                           csp_enum.AZGARD_VALUES_KEY: None,
                           csp_enum.AZGARD_DEBUG: AZGARD_DEBUG,
                           csp_enum.TRANSFORMATION: False}
    return specific_parameters


def _specific_parameters_regression_predictive_model():
    csp_enum = ComponentSpecificParametersEnum()
    transf_type = TransformationTypeEnum()
    descriptor_types = DescriptorTypesEnum()
    specific_parameters = {csp_enum.HIGH: 9,
                           csp_enum.LOW: 4,
                           csp_enum.K: 0.25,
                           csp_enum.TRANSFORMATION: True,
                           csp_enum.TRANSFORMATION_TYPE: transf_type.SIGMOID,
                           csp_enum.SCIKIT: "regression",
                           csp_enum.DESCRIPTOR_TYPE: descriptor_types.ECFP_COUNTS}
    return specific_parameters


def _specific_parameters_classifiaction_predictive_model():
    csp_enum = ComponentSpecificParametersEnum()
    descriptor_types = DescriptorTypesEnum()
    specific_parameters = {csp_enum.HIGH: 9,
                           csp_enum.LOW: 4,
                           csp_enum.K: 0.25,
                           csp_enum.TRANSFORMATION: False,
                           csp_enum.SCIKIT: "classification",
                           csp_enum.DESCRIPTOR_TYPE: descriptor_types.ECFP_COUNTS}
    return specific_parameters


def create_custom_alerts_configuration():
    custom_alerts_list = [
        '[*;r7]',
        '[*;r8]',
        '[*;r9]',
        '[*;r10]',
        '[*;r11]',
        '[*;r12]',
        '[*;r13]',
        '[*;r14]',
        '[*;r15]',
        '[*;r16]',
        '[*;r17]',
        '[#8][#8]',
        '[#6;+]',
        '[#16][#16]',
        '[#7;!n][S;!$(S(=O)=O)]',
        '[#7;!n][#7;!n]',
        'C#C',
        'C(=[O,S])[O,S]',
        '[#7;!n][C;!$(C(=[O,N])[N,O])][#16;!s]',
        '[#7;!n][C;!$(C(=[O,N])[N,O])][#7;!n]',
        '[#7;!n][C;!$(C(=[O,N])[N,O])][#8;!o]',
        '[#8;!o][C;!$(C(=[O,N])[N,O])][#16;!s]',
        '[#8;!o][C;!$(C(=[O,N])[N,O])][#8;!o]',
        '[#16;!s][C;!$(C(=[O,N])[N,O])][#16;!s]']
    sf_enum = ScoringFunctionComponentNameEnum()
    parameters = ComponentParameters(component_type=sf_enum.CUSTOM_ALERTS,
                                        name="custom_alerts",
                                        weight=1.,
                                        smiles=custom_alerts_list,
                                        model_path="",
                                        specific_parameters={})
    return parameters
