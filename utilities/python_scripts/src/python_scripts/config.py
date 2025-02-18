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

            # Set the property on the class if one hasn't already been defined
            try:
                getattr(cls, public_name)
            except AttributeError:
                setattr(cls, public_name, setter)

    return cls


