import io
import json
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np
from typing import List, Tuple

from reinvent_scoring.scoring.utils import _is_development_environment

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.console_invoked.base_console_invoked_component import BaseConsoleInvokedComponent


class AZgard(BaseConsoleInvokedComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self._configuration_path = self.parameters.specific_parameters[self.component_specific_parameters.AZGARD_CONFPATH]
        self._executor_script_path = self.parameters.specific_parameters[self.component_specific_parameters.AZGARD_EXECUTORSCRIPTPATH]
        self._environment_path = self.parameters.specific_parameters[self.component_specific_parameters.AZGARD_ENVPATH]
        self._values_key = self.parameters.specific_parameters[self.component_specific_parameters.AZGARD_VALUES_KEY]

    def _add_debug_mode_if_selected(self, command):
        if self.parameters.specific_parameters.get(self.component_specific_parameters.AZGARD_DEBUG, False)\
                or _is_development_environment():
            command = ' '.join([command, "-debug"])
        return command

    def _create_command(self, step, input_json_path: str, output_json_path: str):
        # use "step" as well for the write-out
        global_variables = "".join(["\"input_json_path:", input_json_path, "\" ",
                                    "\"output_json_path:", output_json_path, "\" ",
                                    "\"step_id:", str(step), "\""])
        command = ' '.join([self._environment_path,
                            self._executor_script_path,
                            "-conf", self._configuration_path,
                            "--global_variables", global_variables])

        # check, if AZgard is to be executed in debug mode, which will cause its loggers to print out
        # much more detailed information
        command = self._add_debug_mode_if_selected(command)
        return command

    def _prepare_input_data_JSON(self, path: str, smiles: List[str]):
        """Needs to look something like:
           {
               "names": ["0", "1", "3"],
               "smiles": ["C#CCCCn1...", "CCCCn1c...", "CC(C)(C)CCC1(c2..."]
           }"""
        names = [str(idx) for idx in range(len(smiles))]
        input_dict = {"names": names,
                      "smiles": smiles}
        with open(path, 'w') as f:
            json.dump(input_dict, f, indent=4)

    def _select_values(self, data: dict) -> list:
        for value_dict in data["results"]:
            if self._values_key == value_dict["value_key"]:
                return value_dict["values"]
        return []

    def _parse_output_data_json(self, path: str) -> Tuple[List[str], List[float]]:
        """Needs to look something like:
           {
               "results": [{
                   "value_key": "docking_score",
                   "values": ["-5.88841", "-5.72676", "-7.30167"]},
                           {
                   "value_key": "shape_similarity",
                   "values": ["0.476677", "0.458017", "0.510676"]},
                           {
                   "value_key": "esp_similarity",
                   "values": ["0.107989", "0.119446", "0.100109"]}],
               "names": ["0", "1", "2"]
           }"""
        names_list = []
        values_list = []

        if not os.path.isfile(path):
            raise FileNotFoundError(f"Output file {path} does not exist, indicating that execution of AZgard failed entirely. Check your setup and the log file.")

        with open(path, 'r') as f:
            data = f.read().replace("\r", "").replace("\n", "")
        data = json.loads(data)
        raw_values_list = self._select_values(data=data)

        for idx in range(len(data["names"])):
            names_list.append(data["names"][idx])
            try:
                score = float(raw_values_list[idx])
            except ValueError:
                score = 0
            values_list.append(score)

        return names_list, values_list

    def _execute_command(self, command: str, final_file_path: str = None):
        # execute the pre-defined command
        with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              shell=True) as proc:
            wrapped_proc_in = io.TextIOWrapper(proc.stdin, "utf-8")
            wrapped_proc_out = io.TextIOWrapper(proc.stdout, "utf-8")
            wrapped_proc_in.close()
            wrapped_proc_out.close()
            proc.wait()
            proc.terminate()

        # wait in case the final file has not been written yet or is empty (in case of a filesystem delay / hick-up)
        if final_file_path is not None:
            for _ in range(5):
                if os.path.isfile(final_file_path) and os.path.getsize(final_file_path) > 0:
                    break
                else:
                    time.sleep(3)

    def _calculate_score(self, smiles: List[str], step) -> np.array:
        # make temporary folder and set input and output paths
        tmp_dir = tempfile.mkdtemp()
        input_json_path = os.path.join(tmp_dir, "input.json")
        output_json_path = os.path.join(tmp_dir, "output.json")

        # save the smiles in an AZgard compatible JSON file
        self._prepare_input_data_JSON(path=input_json_path, smiles=smiles)

        # create the external command
        command = self._create_command(step=step,
                                       input_json_path=input_json_path,
                                       output_json_path=output_json_path)

        # execute the AZgard component
        self._execute_command(command=command, final_file_path=output_json_path)

        # TODO: smiles_ids are not used yet; will be important to match the return values back
        # parse the output
        smiles_ids, scores = self._parse_output_data_json(path=output_json_path)

        # apply transformation
        transformed_scores = self._transformation_function(scores, self.parameters.specific_parameters)

        # clean up
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)

        return np.array(transformed_scores), np.array(scores)
