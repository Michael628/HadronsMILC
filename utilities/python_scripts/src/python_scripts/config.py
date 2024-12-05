import typing as t
from enum import Enum, auto

import nanny.config


class Gamma(Enum):
    ONELINK = auto()
    LOCAL = auto()
    VEC_ONELINK = auto()
    VEC_LOCAL = auto()
    PION_LOCAL = auto()


class ConfigBase:
    pass


def create_config(params: t.Dict) -> ConfigBase:
    """Processes dictionary into object with corresponding key names
    as properties. nested dictionaries are flattened with `_` connecting
    inner and outer keys."""


    def process_val(val: t.Union[str, t.List[str]]):
        """Might use to do some type checking on inputs later"""
        if isinstance(val, t.Dict):
            raise ValueError("No nested dictionaries of depth > 1 in  run config")
        return val

    instance = ConfigBase()
    for k, v in params.items():
        if isinstance(v, t.Dict):
            for k_inner, v_inner in v.items():
                setattr(instance, f"{k}_{k_inner}", process_val(v_inner))
        else:
            setattr(instance, f"{k}", process_val(v))

    return instance


def get_config(config_label: str):
    configs = {
        "epack": nanny.config.create_epack_config,
        "meson": nanny.config.create_op_list_config,
        "high_modes": nanny.config.create_op_list_config,
        'lmi_params': create_config
    }

    if config_label in configs:
        return configs[config_label]
    else:
        raise ValueError(f"No config implementation for `{config_label}`.")
