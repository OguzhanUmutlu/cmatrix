# cmatrix

A simple matrix library made in C++

## Installation

Just install the matrix.cpp file and include it in your file like so:

```cpp
#include "./matrix.cpp"
```

## Creating a matrix

The type syntax for the matrix is `Matrix<elementType, width, height>`

Element type can be any kind of number, or anything, really. It could even be a Vector3 class that implements basic
number operations: + - * /

```cpp
Matrix<double, 3, 3> myMatrix;

std::cout << myMatrix << std::endl;
```

This creates a 3x3 matrix filled with zeros. The output will be:

```cpp
0     0     0
0     0     0
0     0     0
```

## Creating a matrix with a value

```cpp
Matrix<double, 3, 3> myMatrix = {
    3, 5, 8,
    7, 3, 0,
    2, 1, 4
};

std::cout << myMatrix << std::endl;
```

This creates and initializes a matrix with the given values. The output will be:

```cpp
3     5     8
7     3     0
2     1     4
```

## Accessing and setting indices after creation

You can access an index of a matrix using the `myMatrix(i, j)` syntax.

```cpp
Matrix<double, 3, 3> myMatrix = {
    3, 5, 8,
    7, 3, 0,
    2, 1, 4
};

std::cout << myMatrix(0, 1) << std::endl;

myMatrix(0, 1) = 5;

std::cout << myMatrix(0, 1) << std::endl;
```

This will output `7` and `5`.

## Mapping the matrix into another matrix

You can give a function argument to the matrix to map it.

There are two possible kinds of functions you can pass in:

- `myMatrix([](auto x) { return x * x; })`: Squares every value.
- `myMatrix([](auto x, auto i, auto j) { return x + i + j; })`: Increases every value by their positions' sum.

You can optionally give in a second argument as `true` if you want the operations to be applied to the current matrix,
otherwise it will make a replica.

```cpp
Matrix<double, 3, 3> myMatrix = {
    3, 5, 8,
    7, 3, 0,
    2, 1, 4
};

auto mappedMatrix = myMatrix([](double x) { return x * x; });

std::cout << mappedMatrix << std::endl;
```

This will output:

```cpp
9    25    64
49     9     0
4     1    16
```

## Basic operations (=, +, -, *, /, ^, ==, !=)

Here's every use of the operations:

```cpp
Matrix<double, 2, 2> m1 { 1, 2, 3, 4 };
Matrix<double, 2, 2> m2 { 5, 6, 7, 8 };

std::cout << m1 + m2 << std::endl;
m1 += m2;
// 6     8
// 10    12

std::cout << m1 - m2 << std::endl;
m1 -= m2;
// 1     2
// 3     4

std::cout << m1 * m2 << std::endl;
m1 *= m2;
// 23    31
// 34    46

std::cout << m1 / m2 << std::endl;
m1 /= m2;
// 10  -4.5
// 14  -6.5

std::cout << (m1 ^ -1) << std::endl;
// 3.25  -2.25
//    7     -5

std::cout << (m1 ^ 2) << std::endl;
//    37     49
// -15.75  -20.75

std::cout << (m1 == m2) << std::endl; // false
std::cout << (m1 != m2) << std::endl; // true
```

## Converting a matrix to a string

Just cast it.

```cpp
Matrix<double, 2, 2> myMatrix = { 1, 2, 3, 4 };

string matrixStr = (string) myMatrix;
```