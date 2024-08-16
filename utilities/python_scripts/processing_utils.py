import numpy as np
import itertools

def build_format_dict(format_string: str, source_dict: dict):
    """Builds a dictionary from elements of `source_dict` of required variables
    in `format_string`.

    Parameters
    ----------
    format_string : str
        A string with named format variables that need replacement

    source_dict : dict
        The dictionary from which the returned dict is built.

    Returns
    -------
    sub_dict : dict
        Dictionary entries of `source_dict` matching all keyword elements found in `format_string`.
    """

    key_list = [k[1] for k in Formatter().parse(format_string) if k is not None]

    sub_dict = {}
    for key in key_list:
        if key in source_dict and type(source_dict[key]) is str:
            sub_dict[key] = source_dict[key]

    return sub_dict

def build_traversal_list(input_dict: dict):
    """Traverses dictionary values for lists of strings and
    builds a list of all combinations of list values

    Parameters
    ----------
    input_dict : dict
        A dictionary who's values can be typed to str or list[str]

    Returns
    -------
    keys, prod: (list[str], list[tuple[str]])
        `keys` is the list of keys of `input_dict`
        `prod` is a list of all combinations of the values of `input_dict`
    """

    for k,v in input_dict.items():
        if type(v) is not list:
            input_dict[k] = [str(v)]
        else:
            input_dict[k] = [str(v[i]) for i in range(len(v))]

    return input_dict.keys(), list(itertools.product(*input_dict.values()))

def dictToCorrHelper(index, corrOut, corrDict):
    
    # If we've reached the lowest level of the dictionary, copy contents and increment counter
    if type(corrDict) is not dict:
        corrOut[index,:] = corrDict
        index += 1
        return index
    
    # Otherwise recurse over next lower level in dictionary
    for k,v in sorted(corrDict.items()):
        index = dictToCorrHelper(index,corrOut,corrDict[k])
    return index

# Recursively traverses a dictionary,
# creating a multidimensional numpy array
def dictToCorr(corrDict):
    dictShape = tuple()

    # Determine shape of correlator based on number of keys in each nested dictionary
    d = corrDict
    while type(d) is dict:
        dictShape = dictShape + (len(d.keys()),)
        d = list(d.values())[0]
    innerShape = d.shape
    
    dictSize = np.prod(dictShape)
    corrOut = np.empty((dictSize,)+innerShape,dtype=np.complex128)

    # Recursively fill correlator with nested dictionary entries
    index = 0
    for v in corrDict.values():
        index = dictToCorrHelper(index,corrOut,v)
    
    # Reshape to multidimensional structure
    corrOut = np.reshape(corrOut,dictShape+innerShape)
    return corrOut
