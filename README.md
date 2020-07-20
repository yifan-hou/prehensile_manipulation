# Install

# Dependencies
## qhull:
Pull from https://github.com/qhull/qhull
Then simply
```
cd qhull
mkdir build & cd build
cmake ..
make
sudo make install
```

## CDDLIB
Installation:
https://github.com/cddlib/cddlib#build-the-latest-released-version


## CPPLibrary
https://github.com/yifan-hou/cpplibrary

## Pybind11

# Build
Inside the cpp folder,
```
mkdir build
cd build
cmake ..
make
```
A python library is now in build/ folder.

# Run
1. First, launch Matlab, change matlab directory to prehensile_manipulation/matlab.
2. A python test script is located under prehensile_manipulation/python.
```
ipython test2d.py
```

