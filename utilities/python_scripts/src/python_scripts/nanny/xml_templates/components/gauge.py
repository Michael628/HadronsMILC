from python_scripts.nanny.xml_templates.hadrons_templates import \
    module_templates
from python_scripts.utils import deep_copy_dict


def gauge_dp_params(ens: str, series: str,
                    include_unsmeared: bool = False) -> list:

    module_names = ["gauge_fat", "gauge_long"]
    file_prefixes = ['fat', 'lng']

    if include_unsmeared:
        module_names.append("gauge")
        file_prefixes.append('l')

    filenames = [
        "lat/scidac/{smear}{ens}{series}.ildg".format(
            smear=prefix,
            ens=ens,
            series=series
        ) for prefix in file_prefixes
    ]

    modules = [
        deep_copy_dict(module_templates["load_gauge"])
        for n in range(len(module_names))
    ]

    for name, file, module in zip(module_names, filenames, modules):
        module["id"]["name"] = name
        module["options"]["file"] = file

    return modules


def gauge_sp_params(field_names: list[str]) -> list:

    modules = []

    for name in field_names:
        module = deep_copy_dict(module_templates["cast_gauge"])
        module["id"]["name"] = name+"f"
        module["options"]["field"] = name
        modules.append(module)

    return modules


def action_params(name: str, mass: str, fat_label: str,
                  long_label: str, single_precision=False) -> list:

    if single_precision:
        template_name = 'action_float'
    else:
        template_name = 'action'

    module = deep_copy_dict(module_templates[template_name])
    module["id"]["name"] = name
    module["options"]["mass"] = mass
    module["options"]["gaugefat"] = fat_label
    module["options"]["gaugelong"] = long_label

    return [module]


def op_params(name: str, action: str) -> list:
    module = deep_copy_dict(module_templates['op'])
    module["id"]["name"] = name
    module["options"]["action"] = action

    return [module]
