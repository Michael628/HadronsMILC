from python_scripts import ConfigBase
import typing as t


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
            setattr(instance, f"{k}",
                dict(
                    (k_inner,process_val(v_inner))
                    for k_inner, v_inner in v.items()
                )
            )

        elif isinstance(v, t.List):
            setattr(instance, f"{k}", list(map(process_val,v)))
        else:
            setattr(instance, f"{k}", process_val(v))

    return instance
