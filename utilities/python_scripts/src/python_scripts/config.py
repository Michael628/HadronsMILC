import typing as t

T = t.TypeVar('T', bound='ConfigBase')
class ConfigBase:

    @classmethod
    def create(cls: t.Type[T], **kwargs) -> T:
        """Creates a new instance of ConfigBase from a dictionary."""
        return cls(**kwargs)

    @property
    def string_dict(self):
        """Converts all attributes to strings or lists of strings.
        Dictionary attributes are removed from output.
        Returns a dictionary keyed by the attribute labels
        """
        res = {}
        for k, v in self.__dict__.items():
            if k.startswith("_") or isinstance(v, t.Dict):
                continue
            elif isinstance(v,t.List):
                res[k] = list(map(str,v))
            elif isinstance(v,bool):
                res[k] = str(v).lower()
            else:
                if v is not None:
                    res[k] = str(v)

        return res


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
