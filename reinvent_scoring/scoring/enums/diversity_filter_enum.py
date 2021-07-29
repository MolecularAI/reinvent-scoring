class DiversityFilterEnum():
    IDENTICAL_TOPOLOGICAL_SCAFFOLD = "IdenticalTopologicalScaffold"
    IDENTICAL_MURCKO_SCAFFOLD = "IdenticalMurckoScaffold"
    SCAFFOLD_SIMILARITY = "ScaffoldSimilarity"
    NO_FILTER = "NoFilter"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    def __setattr__(self, key, value):
        raise ValueError("Do not assign value to a ScaffoldFilterEnum field.")