from dataclasses import dataclass


@dataclass
class DiversityFilterParameters:
    name: str
    minscore: float
    bucket_size: int
    minsimilarity: float
