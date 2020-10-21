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

# Modeling convention
## Variables
generalized velocity&force:
```
v = [v_HO, v_HH]
f = [f_HO, f_HH]
```
All defined in the hand frame.
v_HO is the object spatial velocity described in the hand frame.
v_HH is the hand body velocity.
f_HO, f_HH are corresponding forces.

## 2D conventions
```
    cone edge 1,      cone edge 2
         -----\-------/------
         |     \  N  /      |
         |      \ ^ /       |
 T_e <---|       \|/        | ---> right sliding
     =============|==================
```
The plane is assumed to be XY plane. So Z vector is perpendicular to the whole plane.
### Mode
* 0 separation
* 1 fixed
* 2 right sliding, left friction force
* 3 left sliding, right friction force

### Jacobian
Normal points inside of the object.
"Left" and "right" are defined locally with respect to the contact normal.
"Left" is computed by cross(z, normal).
```
J = [J_e, 0; -J_h, J_h]:
J_e = [N_e; T_e]
J_h = [N_h; T_h]
```
Each row of N_e, N_h corresponds to a contact normal.
Each row of T_e, T_h corresponds to a contact left direction.

