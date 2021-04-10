

class ScoringFunctionNameEnum:
    CUSTOM_PRODUCT = "custom_product"
    CUSTOM_SUM = "custom_sum"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    def __setattr__(self, key, value):
        raise ValueError("Do not assign value to a ScoringFunctionNameEnum field.")
