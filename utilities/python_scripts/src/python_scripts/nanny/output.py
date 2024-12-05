from dataclasses import dataclass
from python_scripts.config import ConfigBase
import typing as t

@dataclass
class FileOutputConfig(ConfigBase):


def create_epack_output(files: t.Dict):
    pass

def create_meson_output():
    pass

def create_high_modes_output():
    pass


def get_output(output_label: str):
    output = {
        "epack": create_epack_output,
        "meson": create_meson_output,
        "high_modes": create_high_modes_output,
    }

    if output_label in output:
        return output[output_label]
    else:
        raise ValueError(f"No config implementation for `{output_label}`.")
