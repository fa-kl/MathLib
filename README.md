# MathLib
A custom C++ library for mathematics - currently under development. 

## Feature Overview
### Implemented Features
  - Floating point comparison
  - Unit conversions
  - Complex numbers (```Complex```)
  - 2D & 3D Vectors (```Vector2D``` ```Vector3D```)
    - Getters and setters: ```.x()```, ```.y()```, ```.z()```, ```operator[]``` (0-based, C-like), ```operator()``` (1-based, pos. & neg. indices, MATLAB-like) 
    - Element-wise addition, subtraction, multiplicaiton and division (via operators ```+```, ```-```, ```*```)
    - Vector-Vector and vector-scalar addition and subtraction (via operators ```+```, ```-```)
    - Dot- (```*``` and ```dot()```) and Cross-product (```cross()```)
    - Vector-Scalar division (via operator ```/```)
    - Norm (```norm()```), Angles (```angle()```), Normalization (```normalize()```)
    - Exponential and Logarithm (```exp()```, ```log()```)
    - Element-wise square root and pow (```sqrt()```, ```pow()```)
    - Element-wise trigonometric functions (```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```)
    - Conversion to ```std::string``` and printing via ```<iostream>```
    

### Planned Features
  - n-dimensional vectors
  - Matrices
  - More unit conversions
  - Linear algebra, i.e. Eigenvalues, Eigenvectors, SVD, ...

## Build Instructions
### Ubuntu 24.04.3 LTS
Assuming a suitable compiler is installed:
```shell
mkdir build
cd build
cmake ..
make
```

## Documentation
To be done... 
At the moment only doc-strings exist.

## Contribute
This is more a fun project to improve my C++ skills. If you want to improve your C++ skills or have any tipps feel free to contact me! :)
