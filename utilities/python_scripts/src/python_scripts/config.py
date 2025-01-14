import typing as t
from dataclasses import dataclass, fields


def dataclass_with_getters(cls):
    """For any dataclass fields with a single underscore prefix,
    provides a getter property with the underscore removed.

    Note
    ----
    To be used as class decorator."""

    # Apply the dataclass transformation
    cls = dataclass(cls)

    private_fields = [field
                      for field in fields(cls)
                      if field.name.startswith('_') or '__' in field.name]

    # Add properties for each field
    for field in private_fields:
        private_name = field.name

        dunder =  '__' in field.name
        if dunder:
            public_name = field.name.split('__')[1]
        else:
            public_name = field.name.lstrip("_")

        # Add a getter
        getter = property(lambda self, n=private_name: getattr(self, n))

        # Add a setter for properties with only a single underscore
        if not dunder:
            setter = getter.setter(
                lambda self, value, n=private_name: setattr(self, n, value)
            )

            # Set the property on the class
            setattr(cls, public_name, setter)
        else:
            # Set the property on the class
            setattr(cls, public_name, getter)

    return cls


T = t.TypeVar('T', bound='ConfigBase')
class ConfigBase:

    @classmethod
    def create(cls: t.Type[T], **kwargs) -> T:
        """Creates a new instance of ConfigBase from a dictionary.

        Note
        ----
        Checks for dataclass fields stored in the class object.
        If kwargs match a dataclass field apart from a leading underscore,
        e.g. kwargs['mass'] and _mass, the value gets assigned to the underscored object attribute.
        Any kwargs with a double underscore, `dunder_vars`, are presumed to have come from a subclass
        `create` call and are copied directly into the corresponding dunder attribute.
        """

        conflicts = [k for k in kwargs.keys() if k in kwargs and f"_{k}" in kwargs and not k.startswith('_')]

        if any(conflicts):
            raise ValueError(f"Conflict in parameters. Both {conflicts[0]} and _{conflicts[0]} passed to {cls} `create`.")

        class_vars = [f.name for f in fields(cls)]
        obj_vars = {}
        new_vars = {}
        for k,v in kwargs.items():
            if k in class_vars:
                obj_vars[k] = v
            elif f"_{k}" in class_vars and not k.startswith('_'):
                obj_vars[f"_{k}"] = v
            elif k in cls.__dict__:
                raise ValueError(f"Cannot overwrite existing {cls} param, `{k}`. Try relabeling `{k}`")
            else:
                new_vars[k] = v

        obj = cls(**obj_vars)

        for k,v in new_vars.items():
            setattr(obj,k,v)

        return obj

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

    @property
    def public_dict(self):
        """Converts all attributes to strings or lists of strings.
        Dictionary attributes are removed from output.
        Returns a dictionary keyed by the attribute labels
        """
        res = {}
        for k, v in self.__dict__.items():
            if not k.startswith("_"):
                res[k] = v

        return res
