#ifndef DENSE_EIGEN_H
#define DENSE_EIGEN_H

#include <cstddef>
#include <vector>

namespace lobpcg_detail {

class DenseMatrix {
public:
    DenseMatrix() = default;

    DenseMatrix(int rows, int cols) {
        resize(rows, cols);
    }

    void resize(int rows, int cols) {
        rows_ = rows;
        cols_ = cols;
        const std::size_t required = static_cast<std::size_t>(rows_) * cols_;
        if (data_.size() < required) {
            data_.resize(required, 0.0);
        }
    }

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    double& operator()(int r, int c) {
        return data_[static_cast<std::size_t>(r) * cols_ + c];
    }

    double operator()(int r, int c) const {
        return data_[static_cast<std::size_t>(r) * cols_ + c];
    }

    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }

private:
    int rows_ = 0;
    int cols_ = 0;
    std::vector<double> data_;
};

class SmallGeneralizedEigenSolver {
public:
    virtual ~SmallGeneralizedEigenSolver() = default;

    virtual void solve(const DenseMatrix& S, const DenseMatrix& T,
                       std::vector<double>& eigenvalues,
                       DenseMatrix& eigenvectors) const = 0;
};

class JacobiGeneralizedEigenSolver final : public SmallGeneralizedEigenSolver {
public:
    void solve(const DenseMatrix& S, const DenseMatrix& T,
               std::vector<double>& eigenvalues,
               DenseMatrix& eigenvectors) const override;
};

} // namespace lobpcg_detail

#endif // DENSE_EIGEN_H
