from dataclasses import dataclass
from typing import List


@dataclass
class ScoringFuncionParameters:
    name: str
    parameters: List[dict]
    parallel: bool = False
