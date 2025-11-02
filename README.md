# MathLib
A custom C++ library for mathematics - currently under development. 

## Features
  - Floating-point comparison
  - Unit conversions
  - Complex numbers
  - n-dimensional vectors
  - Matrices of arbitrary size
    
### Floating-Point Comparison
  - Comparison for equality: ```isFuzzyEqual()```
  - Check for strictly greater: ```isStrictFuzzyGreater()```
  - Check for greater or equal: ```isFuzzyGreater()```
  - Check for strictly smaller: ```isStrictFuzzySmaller()```
  - Check for smaller or equal: ```isFuzzySmaller()```
    
### Unit Conversions
  - Radiants to degrees: ```rad2deg()```
  - Degrees to radiants: ```deg2rad()```
    
### Complex Numbers (```class Complex```)
  - Getters: ```real()``` and ```imag()```
  - Arithmetic operations: ```+```, ```-```, ```*```, ```/```, ```+=```, ```-=```, ```*=```, ```/=```
    - ```Complex``` & ```Complex```
    - ```Complex``` & ```real_t```
    - ```real_t``` & ```Complex```
  - Magnitude: ```abs()```
  - Squared magnitude: ```abs2()```
  - Argument: ```arg()```
  - Conjungate: ```conj()```
  - Exponential map: ```exp()```
  - Logarithmic map: ```log()```
  - Trigonometric maps: ```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```
  - Power: ```pow()```
  - Square root: ```sqrt()```
  - Conversion to ```std::string```: ```to_string()```
  - Serialization and printing to the default output stream: ```<<``` and ```print()```

### n-Dimensional Vectors (```Vector```)
  - Constructors:
    - Empty vector: ```Vector()```
    - n-dimensional vector: ```Vector(n)```
    - Copy & move constructors and operators
    - Vector of all ones: ```ones(n)```
    - Vector of all zeros: ```zeros(n)```
  - Element access
    - 0-based: ```operator[]```
    - 1-based: ```operator()``` allows positive & negative indices
  - Information acess
    - Number of elements: ```length()```
    - Dimension: ```size()```
    - Number of rows: ```rows()```
    - Number of columns: ```cols()```
  - Element-wise addition and subtraction
    - ```Vector``` & ```Vector```: ```+```, ```-```
    - ```Vector``` & scalar: ```+```, ```-```
    - scalar & ```Vector```: ```+```, ```-```
  - Element-wise multiplicaiton and division
    - Vector-Vector: ```elmul()```, ```eldiv()```
    - Vector-Scalar: ```*```, ```/```
    - Scalar-Vector: ```*```, ```/```
  - Dot product: ```*```, ```dot()```
  - Cross product: ```cross()```
  - Outer product: ```outer()```
  - Transpose: ```transpose()```
  - Norm: ```norm()```
  - Normalization: ```normalize()```
  - Sum of elements: ```sum()```
  - Difference between vector elements: ```diff()```
  - Element-wise rounding: ```round()```
  - Element-wise absolute values: ```abs()```
  - Maximum and minimum: ```max()```, ```min()```
  - Mean, standard deviation, variance: ```mean()```, ```std()```, ```var()```
  - Exponential and logarithm: ```exp()```, ```log()```
  - Element-wise square root and pow: ```sqrt()```, ```pow()```
  - Element-wise trigonometric functions: ```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```
  - Element-wise logical operators: ```&&```, ```||```, ```^```
  - Element-wise comparison via operators: ```==```, ```!=```, ```>```, ```>=```, ```<```, ```<=```
  - Element-wise floating-point comparison via methods: ```isFuzzyEqual()```, ```isStrictFuzzyGreater()```, ```isFuzzyGreater()```, ```isStrictFuzzySmaller()```, ```isFuzzySmaller()```
  - Conversion to ```std::string``` and printing via ```<iostream>```
  - and more ...
    
### Matrices (```Matrix```)  
  - Constructors:
    - Empty matrix: ```Matrix()```
    - (m, n)-dimensional matrix: ```Matrix(m, n)```
    - Copy & move constructors and operators
    - Identity matrix: ```eye(n)```
    - Matrix of all ones: ```ones(m, n)```
    - Matrix of all zeros: ```zeros(m, n)```
    - Diagonal matrix: ```diag(Vector)``` or ```diag({a, b, ...})```
  - Element access
    - 0-based: ```operator[]``` (linear access only)
    - 1-based: ```operator()``` (positive & negative linear and positive & negative 2D indices)
  - Information acess
    - Number of elements: ```length()```
    - Dimension: ```size()```
    - Number of rows: ```rows()```
    - Number of columns: ```cols()```
    - Rank: ```rank()```
  - Element-wise addition and subtraction
    - ```Matrix``` & ```Matrix```: ```+```, ```-```
    - ```Matrix``` & scalar: ```+```, ```-```
    - scalar & ```Matrix```: ```+```, ```-```
  - Element-wise multiplicaiton
    - Matrix-Matrix: ```*```
    - Matrix-Vector: ```*```
    - Vector-Matrix: ```*```
    - Matrix-Scalar: ```*```
    - Scalar-Matrix: ```*```
  - Element-wise division
    - Matrix-scalar: ```/```
  - Dot product: ```*```, ```dot()```, ```colDot()```, ```rowDot()```
  - Outer product: ```outer()```
  - Transpose: ```transpose()```
  - Trace:```trace()``` 
  - Matrix norms: ```norm()```, ```vecnorm()```
  - Orthogonalization: ```orth()```
  - QR decomposition: ```qr()```
  - Eigenvalues and -vectors: ```eig()```
  - SVD: ```svd()```
  - Normalization: ```normalize()``` (rows or columns)
  - Sum of elements: ```sum()``` (all, rows-wise, column-wise)
  - Element-wise rounding: ```round()```
  - Element-wise absolute values: ```abs()```
  - Maximum and minimum: ```max()```, ```min()``` (all, rows-wise, column-wise)
  - Mean, standard deviation, variance: ```mean()```, ```std()```, ```var()```  (all, rows-wise, column-wise)
  - Exponential and logarithm: ```exp()```, ```log()```
  - Element-wise square root and pow: ```sqrt()```, ```elpow()```
  - Element-wise trigonometric functions: ```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```
  - Element-wise logical operators: ```&&```, ```||```, ```^```
  - Element-wise comparison via operators: ```==```, ```!=```, ```>```, ```>=```, ```<```, ```<=```
  - Element-wise floating-point comparison via methods: ```isFuzzyEqual()```, ```isStrictFuzzyGreater()```, ```isFuzzyGreater()```, ```isStrictFuzzySmaller()```, ```isFuzzySmaller()```
  - Conversion to ```std::string``` and printing via ```<iostream>```
  - and more ...

## Build Instructions
### Ubuntu 24.04.3 LTS
Assuming a suitable compiler is installed:
```shell
mkdir build
cd build
cmake ..
make
```
or using the provided shell scripts, i.e. ```build.sh``` and ```run.sh```.

## Documentation
To be done... 
At the moment only doc-strings exist.

## Contribute
This is more a fun project to improve my C++ skills. If you want to improve your C++ skills or have any tipps feel free to contact me! :)
