import pickle

from typing import List

from reinvent_scoring.scoring.predictive_model.model_container import ModelContainer
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory
from reinvent_scoring.scoring.enums import TransformationTypeEnum


class PredictivePropertyComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.activity_model = self._load_model(parameters)
        self.transformation = self._assign_transformation(parameters.specific_parameters)

    def calculate_score(self, molecules: List) -> ComponentSummary:
        score, raw_score = self._predict_and_transform(molecules)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters, raw_score=raw_score)
        return score_summary

    def _predict_and_transform(self, molecules: List):
        score = self.activity_model.predict(molecules, self.parameters.specific_parameters)
        transformed_score = self._apply_transformation(score, self.parameters.specific_parameters)
        return transformed_score, score

    def _load_model(self, parameters: ComponentParameters):
        try:
            activity_model = self._load_container(parameters)
        except:
            raise Exception(f"The loaded file {parameters.model_path} isn't a valid scikit-learn model")
        return activity_model

    def _load_container(self, parameters: ComponentParameters):
        with open(parameters.model_path, "rb") as f:
            scikit_model = pickle.load(f)
            packaged_model = ModelContainer(scikit_model, parameters.specific_parameters)
        return packaged_model

    def _apply_transformation(self, predicted_activity, parameters: dict):
        if parameters.get(self.component_specific_parameters.TRANSFORMATION, False):
            activity = self.transformation(predicted_activity, parameters)
        else:
            activity = predicted_activity
        return activity

    def _assign_transformation(self, specific_parameters: dict):
        transformation_type = TransformationTypeEnum()
        """classification models should not have any prediction transformations"""
        if specific_parameters[self.component_specific_parameters.SCIKIT] == "classification":
            specific_parameters[self.component_specific_parameters.TRANSFORMATION] = False
            specific_parameters[self.component_specific_parameters.TRANSFORMATION_TYPE] = transformation_type.NO_TRANSFORMATION
        factory = TransformationFactory()
        transform_function = factory.get_transformation_function(specific_parameters)
        return transform_function
