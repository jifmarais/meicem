#pragma once
#include <armadillo>
#include <complex>
#include <vector>

class ComplexMatrix
{
    public:
        typedef std::complex<double> complex_type;
        typedef std::vector<complex_type> vector_type;
        typedef std::vector<vector_type> matrix_type;

        ComplexMatrix(unsigned rows, unsigned columns, complex_type initialValue);
        ComplexMatrix(unsigned rows, unsigned columns);
        ComplexMatrix(unsigned rows);
        virtual               ~ComplexMatrix();

        complex_type&          operator()(unsigned row, unsigned col);
        const complex_type&    operator()(unsigned row, unsigned col) const;

        static ComplexMatrix   identity(unsigned rows);
        unsigned               getRowCount() const;
        unsigned               getColumnCount() const;
        ComplexMatrix&         operator=(const ComplexMatrix& rhs);
        bool                   operator==(const ComplexMatrix& rhs) const;
        bool                   operator!=(const ComplexMatrix& rhs) const;
        ComplexMatrix          operator+(const ComplexMatrix& rhs) const;
        ComplexMatrix&         operator+=(const ComplexMatrix& rhs);
        ComplexMatrix          operator+(const complex_type& rhs) const;
        ComplexMatrix          operator-(const ComplexMatrix& rhs) const;
        ComplexMatrix&         operator-=(const ComplexMatrix& rhs);
        ComplexMatrix          operator-(const complex_type& rhs) const;
        ComplexMatrix          hadamard(const ComplexMatrix& rhs) const;
        ComplexMatrix          operator*(const ComplexMatrix& rhs) const;
        ComplexMatrix&         operator*=(const ComplexMatrix& rhs);
        ComplexMatrix          operator*(const complex_type& rhs) const;
        ComplexMatrix          operator/(const ComplexMatrix& rhs) const;
        ComplexMatrix&         operator/=(const ComplexMatrix& rhs);
        ComplexMatrix          operator/(const complex_type& rhs) const;
        ComplexMatrix 	       conjugateTranspose() const;
        ComplexMatrix 	       transpose() const;
        ComplexMatrix          subMatrix(unsigned rowStart, unsigned rowEnd, unsigned columnStart, unsigned columnEnd) const;
        ComplexMatrix          minorMatrix(unsigned row, unsigned col) const;
        complex_type           determinant() const;
        ComplexMatrix          inverse();
        void                   print() const;

    protected:
    private:
        void                   swapValues(unsigned rowA, unsigned columnA, unsigned rowB, unsigned columnB);
        void                   init(unsigned rows, unsigned columns, const complex_type& initialValue);

        matrix_type            m_matrix;
        unsigned               m_rows;
        unsigned               m_columns;
};

ComplexMatrix operator+(ComplexMatrix::complex_type d, const ComplexMatrix& rhs);
ComplexMatrix operator-(ComplexMatrix::complex_type d, const ComplexMatrix& rhs);
ComplexMatrix operator*(ComplexMatrix::complex_type d, const ComplexMatrix& rhs);
ComplexMatrix operator/(ComplexMatrix::complex_type d, const ComplexMatrix& rhs);
