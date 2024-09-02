#! /bin/bash

moduleURL=https://files.pythonhosted.org/packages/72/04/331c8d0873fab15309e57b803aa0de23ba074ba16e0239445dbb6c63dfa1/dict2xml-1.7.1.tar.gz

folder=dict2xml-1.7.1
tarFile=dict2xml-1.7.1.tar.gz

wget ${moduleURL}
tar -xzf ${tarFile}

pushd ${folder}
python3 setup.py install --user
popd

rm -rf ${folder}
rm ${tarFile}
