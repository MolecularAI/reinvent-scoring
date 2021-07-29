from dataclasses import dataclass
from typing import List


# TODO: Remove this later
@dataclass
class ScoringFuncionParameters:
    name: str
    parameters: List[dict]
    parallel: bool = False

@dataclass
class ScoringFunctionParameters:
    name: str
    parameters: List[dict]
    parallel: bool = False