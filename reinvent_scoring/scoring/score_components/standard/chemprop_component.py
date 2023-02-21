from contextlib import redirect_stdout
import io

from typing import List

import chemprop
from rdkit import Chem

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_transformations import TransformationFactory
from reinvent_scoring.scoring.enums import TransformationTypeEnum, TransformationParametersEnum


class ChemPropComponent(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        args = [
            '--checkpoint_dir', parameters.specific_parameters['checkpoint_dir'],  # ChemProp models directory
            '--test_path', '/dev/null',  # required
            '--preds_path', '/dev/null'  # required
        ]

        rdkit2d = parameters.specific_parameters.get('rdkit_2d_normalized', False)

        if rdkit2d:
            args.extend([
                '--features_generator', 'rdkit_2d_normalized',  '--no_features_scaling'
            ])

        f = io.StringIO()
        with redirect_stdout(f):
            self.chemprop_args = chemprop.args.PredictArgs().parse_args(args)
            self.chemprop_model = chemprop.train.load_model(args=self.chemprop_args)

        self._transformation_function = self._assign_transformation(parameters.specific_parameters)

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        smilies = [[Chem.MolToSmiles(molecule)] if molecule else ['INVALID'] for molecule in molecules]

        preds = chemprop.train.make_predictions(
            model_objects=self.chemprop_model,
            smiles=smilies,
            args=self.chemprop_args,
            return_invalid_smiles=True,
            return_uncertainty=False
        )
        raw_scores = [val[0] if not 'Invalid SMILES' in val else 0.0 for val in preds]  # FIXME: score 0.0
        scores =  self._apply_transformation(raw_scores, self.parameters.specific_parameters)

        score_summary = ComponentSummary(total_score=scores, parameters=self.parameters, raw_score=raw_scores)

        return score_summary

    def _apply_transformation(self, predicted_activity, parameters: dict):
        transform_params = parameters.get(self.component_specific_parameters.TRANSFORMATION)

        if transform_params:
            activity = self._transformation_function(predicted_activity, transform_params)
        else:
            activity = predicted_activity

        return activity

    def _assign_transformation(self, specific_parameters: dict):
        transformation_type = TransformationTypeEnum()
        transform_params = specific_parameters.get(self.component_specific_parameters.TRANSFORMATION)

        if not transform_params:
            specific_parameters[self.component_specific_parameters.TRANSFORMATION] = {
                    TransformationParametersEnum.TRANSFORMATION_TYPE: transformation_type.NO_TRANSFORMATION
                }

        factory = TransformationFactory()

        return factory.get_transformation_function(transform_params)



if __name__ == '__main__':
    from rdkit import Chem

    cp = ComponentParameters('', '', 1.0, {
            'checkpoint_dir': '3CLPro_6w63',
            'rdkit_2d_normalized': True,
            'transformation': {
                'transformation_type': 'reverse_sigmoid',
                'high': -5.0,
                'low': -35.0,
                'k': 0.4
            }
        })

    component = ChemPropComponent(cp)

    smilies = [
        'CC(C)(C)NC(=O)[C@H]1N(C(=O)[C@@H](O)[C@H](Cc2ccccc2)NC(=O)c2ccc(Cl)cc2)CSC1(C)C',  # CHEMBL329432
        'INVALID',
        'Cc1nc2n(c(=O)c1CCN1CCC(c3noc4cc(F)ccc34)CC1)CCCC2'  # CHEMBL85, RISPERIDONE
    ]
    mols = [Chem.MolFromSmiles(smiles) for smiles in smilies]
    score_summary = component.calculate_score(mols)
    print(score_summary.total_score)
    print(score_summary.raw_score)
