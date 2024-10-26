#include <iostream>
#include "matrix.cpp"

int main() {
    Matrix<double, 2, 2> m1{1, 2, 3, 4};
    Matrix<double, 2, 2> m2{5, 6, 7, 8};

    std::cout << m1 + m2 << std::endl;
    m1 += m2;

    std::cout << m1 - m2 << std::endl;
    m1 -= m2;

    std::cout << m1 * m2 << std::endl;
    m1 *= m2;

    std::cout << m1 / m2 << std::endl;
    m1 /= m2;

    std::cout << (m1 == m2) << std::endl;
    std::cout << (m1 != m2) << std::endl;

    std::cout << (m1 ^ -1) << std::endl;
    std::cout << (string) (m1 ^ 2) << std::endl;
    string matrixStr = (string) m1;

    return 0;
}
