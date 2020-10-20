# Install

# Dependencies

## PPL
First, install GMP and glpk:
```
sudo apt install libgmp-dev
sudo apt install libglpk-dev
```

Download from https://www.bugseng.com/ppl-download
Extract somewhere,
```
cd ppl-1.2
./configure
make
sudo make install
```
I got some errors when calling "sudo make install", but that turned out to be fine.
To see usage of PPL, look at ppl-1.2/demos/ppl_lcdd and ppl-1.2/tests/Polyhedron

## qhull:
Pull from https://github.com/qhull/qhull
Then simply
```
cd qhull/build
cmake ..
make
sudo make install
```
## CPPLibrary
https://github.com/yifan-hou/cpplibrary

## Pybind11
```
git clone git@github.com:pybind/pybind11.git
cd pybind11
mkdir build
cd build
cmake ..
make check -j 4
sudo make install
```

## Matlab dependencies
### Install Matlab Engine for Python
https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

### Contact mode enumeration
```
git clone git@github.com:XianyiCheng/contact_mode_enumeration_2d.git
```
add to path

### matlablibrary
```
git clone git@github.com:yifan-hou/matlablibrary.git
```
add to path

## Setup file
```shell
source /path/to/prehensile_manipulation/setup.bash
echo 'source /path/to/prehensile_manipulation/setup.bash' >> ~/.bashrc
```

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
2. Run the test script in a terminal. The test script is located under prehensile_manipulation/python.
```
python3 test2d.py
```

