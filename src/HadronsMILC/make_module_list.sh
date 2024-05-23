#!/usr/bin/env bash

echo 'modules_cpp =\' > modules.inc
find Modules -name '*.cpp' -type f -print | LC_ALL=C sort | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
echo '' >> modules.inc
echo 'modules_hpp =\' >> modules.inc
find Modules -name '*.hpp' -type f -print | LC_ALL=C sort | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
echo '' >> modules.inc
rm -f Modules.hpp
rm -f Modules.cpp
for f in $(find Modules -name '*.hpp' | LC_ALL=C sort); do
	echo "#include \"${f}\"" >> Modules.hpp
done

for f in $(find Modules -name '*.cpp' | LC_ALL=C sort); do
	echo "#include \"${f}\"" >> Modules.cpp
done
