#ifndef CMATRIX_MATRIX_CPP
#define CMATRIX_MATRIX_CPP

#include <stdexcept>
#include <functional>
#include <sstream>
#include <iomanip>

using namespace std;

template<typename T, size_t w, size_t h>
class Matrix {
public:
    Matrix() {
        for (int i = 0; i < w * h; ++i) {
            m[i] = T();
        }
    };

    explicit Matrix(const T l[w * h]) {
        for (int i = 0; i < w * h; ++i) {
            m[i] = l[i];
        }
    };

    Matrix(initializer_list<T> values) {
        if (values.size() != w * h) {
            throw invalid_argument("matrix initializer was expected to have exactly " + to_string(w * h) + " elements");
        }

        int i = 0;
        for (const T &value: values) {
            m[i++] = value;
        }
    };

    Matrix(initializer_list<initializer_list<T>> values) {
        if (values.size() != w) {
            throw invalid_argument("matrix initializer was expected to have exactly " + to_string(w) + " elements");
        }

        int i = 0;
        for (const auto &inner_list: values) {
            if (inner_list.size() != h) {
                throw invalid_argument(
                    "matrix initializer's every row was expected to have exactly " + to_string(h) + " elements");
            }

            for (const T &value: inner_list) {
                m[i++] = value;
            }
        }
    };

    Matrix(const Matrix<T, w, h> &other) {
        for (int i = 0; i < w * h; ++i) {
            m[i] = other.m[i];
        }
    };

    T m[w * h];

    // m(x, y) = v;
    T &operator()(size_t x, size_t y) {
        if (x >= w || y >= h) {
            throw out_of_range("index out of range");
        }

        return m[y * w + x];
    };

    // m(x, y)
    T operator()(size_t x, size_t y) const {
        if (x >= w || y >= h) {
            throw out_of_range("index out of range");
        }

        return m[y * w + x];
    };

    // m([](T x) { return sqrt(x); })
    Matrix<T, w, h> operator()(const function<T(T)> &map, bool apply = false) {
        Matrix<T, w, h> res = apply ? *this : Matrix<T, w, h>(*this);
        for (int i = 0; i < w * h; ++i) {
            res.m[i] = map(res.m[i]);
        }

        return res;
    };

    // m([](T x, size_t i, size_t j) { return sqrt(x); })
    Matrix<T, w, h> operator()(const function<T(T, size_t, size_t)> &map, bool apply = false) {
        Matrix<T, w, h> res = apply ? *this : Matrix<T, w, h>(*this);
        for (size_t i = 0; i < w; ++i) {
            for (size_t j = 0; j < h; ++j) {
                size_t k = i + j * w;
                res.m[k] = map(res.m[k]);
            }
        }

        return res;
    };

    // m1 = m2;
    Matrix<T, w, h> &operator=(const Matrix<T, w, h> &other) {
        if (this == &other) return *this;
        for (int i = 0; i < w * h; ++i) {
            m[i] = other.m[i];
        }

        return *this;
    };

    // m1 += m2;
    Matrix<T, w, h> &operator+=(const Matrix<T, w, h> &other) {
        for (int i = 0; i < w * h; ++i) {
            m[i] += other.m[i];
        }

        return *this;
    };

    // m1 -= m2;
    Matrix<T, w, h> &operator-=(const Matrix<T, w, h> &other) {
        for (int i = 0; i < w * h; ++i) {
            m[i] -= other.m[i];
        }

        return *this;
    };

    // m1 *= m2;
    Matrix<T, w, h> &operator*=(const Matrix<T, w, h> &other) {
        return *this = *this * other;
    };

    // m1 *= m2;
    Matrix<T, w, h> &operator/=(const Matrix<T, w, h> &other) {
        return *this = *this / other;
    };

    // m *= number;
    Matrix<T, w, h> &operator*=(T number) {
        for (int i = 0; i < w * h; ++i) {
            m[i] *= number;
        }

        return *this;
    };

    // m /= number;
    Matrix<T, w, h> &operator/=(T number) {
        if (number == 0) {
            throw logic_error("division by zero");
        }

        for (int i = 0; i < w * h; ++i) {
            m[i] /= number;
        }

        return *this;
    };

    // m1 + m2
    Matrix<T, w, h> operator+(const Matrix<T, w, h> &other) const {
        Matrix<T, w, h> result;
        for (int i = 0; i < w * h; ++i) {
            result.m[i] = m[i] + other.m[i];
        }

        return result;
    };

    // m1 - m2
    Matrix<T, w, h> operator-(const Matrix<T, w, h> &other) const {
        Matrix<T, w, h> result;
        for (int i = 0; i < w * h; ++i) {
            result.m[i] = m[i] - other.m[i];
        }

        return result;
    };

    // m1 * m2
    template<size_t z>
    Matrix<T, w, z> operator*(const Matrix<T, h, z> &other) {
        Matrix<T, w, z> result;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < z; ++j) {
                result.m[i * z + j] = 0;
                for (int k = 0; k < h; ++k) {
                    result.m[i * z + j] += (*this)(i, k) * other(k, j);
                }
            }
        }

        return result;
    };

    // m1 / m2
    Matrix<T, w, h> operator/(const Matrix<T, w, h> &other) {
        if (w != h) {
            throw logic_error("cannot apply division on non-square matrices");
        }

        return this->operator*(other.inverse());
    }

    // m * number
    Matrix<T, w, h> operator*(T number) const {
        Matrix<T, w, h> result;
        for (int i = 0; i < w * h; ++i) {
            result.m[i] = m[i] * number;
        }

        return result;
    };

    // m / number
    Matrix<T, w, h> operator/(T number) const {
        if (number == 0) {
            throw logic_error("division by zero");
        }

        Matrix<T, w, h> result;
        for (int i = 0; i < w * h; ++i) {
            result.m[i] = m[i] / number;
        }

        return result;
    };

    // m1 == m2
    bool operator==(const Matrix<T, w, h> &other) const {
        for (int i = 0; i < w * h; ++i) {
            if (m[i] != other.m[i]) {
                return false;
            }
        }

        return true;
    };

    // m1 != m2
    bool operator!=(const Matrix<T, w, h> &other) const {
        return !(*this == other);
    };

    // m[i] = number;
    T &operator[](size_t index) {
        if (index >= w * h) {
            throw out_of_range("index out of range");
        }
        return m[index];
    };

    // m[i]
    T operator[](size_t index) const {
        if (index >= w * h) {
            throw out_of_range("index out of range");
        }
        return m[index];
    };

    [[maybe_unused]] [[nodiscard]] pair<size_t, size_t> size() const {
        return {w, h};
    };

    Matrix<T, w, h> operator--() {
        if (w != h) {
            throw logic_error("w != h");
        }

        for (int i = 0; i < w; ++i) {
            for (int j = i + 1; j < w; ++j) {
                swap(m[i * w + j], m[j * w + i]);
            }
        }

        return *this;
    };

    Matrix<T, h, w> transpose() {
        Matrix<T, h, w> transposed;
        for (int row = 0; row < w; ++row) {
            for (int col = 0; col < h; ++col) {
                transposed(col, row) = (*this)(row, col);
            }
        }
        return transposed;
    };

    // ~m
    inline Matrix<T, w, h> operator~() {
        return transpose();
    };

    Matrix<T, w, h> operator^(int x) {
        if (x < 0) return inverse() ^ (-x);
        if (x == 0) return identity();
        if (x == 1) return *this;
        if (x == 2) return *this * *this;

        Matrix<T, w, h> result = *this;
        for (int i = 0; i < x - 1; ++i) {
            result *= *this;
        }

        return result;
    };

    Matrix<T, w, h> &operator^=(int x) {
        return *this = *this ^ x;
    };

    Matrix<T, w, h> inverse() const {
        auto inv = identity();
        Matrix<T, w, h> augmented(*this);

        // https://en.wikipedia.org/wiki/Gaussian_elimination
        for (size_t i = 0; i < w; ++i) {
            T pivot = augmented(i, i);
            if (abs(pivot) < 1e-9) {
                throw runtime_error("singular matrices don't have an inverse");
            }

            for (size_t j = 0; j < w; ++j) {
                augmented(i, j) /= pivot;
                inv(i, j) /= pivot;
            }

            for (size_t k = 0; k < w; ++k) {
                if (k == i) continue;
                T factor = augmented(k, i);
                for (size_t j = 0; j < w; ++j) {
                    augmented(k, j) -= factor * augmented(i, j);
                    inv(k, j) -= factor * inv(i, j);
                }
            }
        }

        return inv;
    };

    explicit operator string() const {
        ostringstream oss;
        oss << *this;
        return oss.str();
    };

    friend ostream &operator<<(ostream &os, const Matrix<T, w, h> &matrix) {
        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                os << setw(5) << matrix(i, j) << " ";
            }
            os << "\n";
        }
        return os;
    };

    static Matrix<T, w, h> identity() {
        Matrix<T, w, h> result;
        for (int i = 0; i < w * h; ++i) {
            result.m[i] = T();
        }
        for (int i = 0; i < min(w, h); ++i) {
            result(i, i) = T(1);
        }
        return result;
    };
};

#endif //CMATRIX_MATRIX_CPP
