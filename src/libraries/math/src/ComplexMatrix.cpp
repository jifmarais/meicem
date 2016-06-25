#include "ComplexMatrix.hpp"
#include <iostream>
#include <math.h>
#include <assert.h>

ComplexMatrix::ComplexMatrix(unsigned rows, unsigned columns, const complex_type initialValue)
{
    m_rows = 0;
    m_columns = 0;
    init(rows, columns, initialValue);
}

ComplexMatrix::ComplexMatrix(unsigned rows, unsigned columns)
{
    m_rows = 0;
    m_columns = 0;
    complex_type defaultValue {0.0, 0.0};
    init(rows, columns, defaultValue);
}

// Create a square matrix
ComplexMatrix::ComplexMatrix(unsigned rows)
{
    m_rows = 0;
    m_columns = 0;
    complex_type defaultValue {0.0, 0.0};
    init(rows, rows, defaultValue);
}

ComplexMatrix::~ComplexMatrix()
{
    //dtor
}

ComplexMatrix ComplexMatrix::identity(unsigned rows)
{
    complex_type zero {0.0, 0.0};
    complex_type one  {1.0, 0.0};
    ComplexMatrix m {rows, rows, zero};

    for (unsigned rr=0; rr < rows; ++rr)
    {
        m(rr, rr) = one;
    }

    return m;
}

void ComplexMatrix::init(unsigned rows, unsigned columns, const complex_type& initialValue)
{
    // If the matrix has been used in the base, clear it completely
    if ( (m_rows > 0) || (m_columns > 0) )
    {
        for (unsigned rr=0; rr < m_rows; ++rr)
        {
            m_matrix[rr].clear();
            m_matrix[rr].shrink_to_fit();
        }
        m_matrix.clear();
        m_matrix.shrink_to_fit();
    }

    // Set the new values
    m_rows = rows;
    m_columns = columns;

    vector_type tmp;
    tmp.resize(columns, initialValue);
    m_matrix.resize(rows, tmp);
}

ComplexMatrix::complex_type& ComplexMatrix::operator()(unsigned row, unsigned col)
{
    assert(row < m_rows);
    assert(col < m_columns);
    return this->m_matrix[row][col];
}

const ComplexMatrix::complex_type& ComplexMatrix::operator()(unsigned row, unsigned col) const
{
    assert(row < m_rows);
    assert(col < m_columns);
    return this->m_matrix[row][col];
}

bool ComplexMatrix::operator==(const ComplexMatrix& rhs) const
{
    if ( &rhs == this )
    {
      return true;
    }

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            if ( m_matrix[rr][cc] != rhs(rr, cc) )
            {
                return false;
            }
        }
    }

    return true;
}

bool ComplexMatrix::operator!=(const ComplexMatrix& rhs) const
{
    return not ( this->operator==(rhs) );
}

ComplexMatrix& ComplexMatrix::operator=(const ComplexMatrix& rhs)
{
    if ( &rhs == this )
    {
      return *this;
    }

    this->init(rhs.getRowCount(), rhs.getColumnCount(), {0.0, 0.0});

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m_matrix[rr][cc] = rhs(rr, cc);
        }
    }

    return *this;
}

ComplexMatrix ComplexMatrix::operator+(const ComplexMatrix& rhs) const
{
    assert(m_rows == rhs.getRowCount());
    assert(m_columns == rhs.getColumnCount());

    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] + rhs(rr, cc);
        }
    }
    return m;
}

ComplexMatrix& ComplexMatrix::operator+=(const ComplexMatrix& rhs)
{
    *this = *this + rhs;
    return *this;
}

ComplexMatrix ComplexMatrix::operator+(const complex_type& rhs) const
{
    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] + rhs;
        }
    }
    return m;
}

ComplexMatrix operator+(ComplexMatrix::complex_type c, const ComplexMatrix& rhs)
{
    return rhs + c;
}

ComplexMatrix ComplexMatrix::operator-(const ComplexMatrix& rhs) const
{
    assert(m_rows == rhs.getRowCount());
    assert(m_columns == rhs.getColumnCount());

    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] - rhs(rr, cc);
        }
    }
    return m;
}

ComplexMatrix& ComplexMatrix::operator-=(const ComplexMatrix& rhs)
{
    *this = *this - rhs;
    return *this;
}

ComplexMatrix ComplexMatrix::operator-(const complex_type& rhs) const
{
    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] - rhs;
        }
    }
    return m;
}

ComplexMatrix operator-(ComplexMatrix::complex_type c, const ComplexMatrix& rhs)
{
    return rhs - c;
}

ComplexMatrix ComplexMatrix::operator*(const ComplexMatrix& rhs) const
{
    assert(m_columns == rhs.getRowCount());

    ComplexMatrix m {m_rows, rhs.getColumnCount()};

    const unsigned rows = m.getRowCount();
    const unsigned cols = m.getColumnCount();
    const unsigned rhsRows = rhs.getRowCount();

    for (unsigned mrr=0; mrr < rows; ++mrr)
    {
        for (unsigned rhsrr=0; rhsrr < rhsRows; ++rhsrr)
        {
            for (unsigned mcc=0; mcc < cols; ++mcc)
            {
                m(mrr,mcc) += this->m_matrix[mrr][rhsrr] * rhs(rhsrr,mcc);
            }
        }
    }

    return m;
}

// Element wise multiplication for matrix
ComplexMatrix ComplexMatrix::hadamard(const ComplexMatrix& rhs) const
{
    assert(m_rows == rhs.getRowCount());
    assert(m_columns == rhs.getColumnCount());

    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] * rhs(rr, cc);
        }
    }
    return m;
}

ComplexMatrix& ComplexMatrix::operator*=(const ComplexMatrix& rhs)
{
    *this = *this * rhs;
    return *this;
}

ComplexMatrix ComplexMatrix::operator*(const complex_type& rhs) const
{
    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] * rhs;
        }
    }
    return m;
}

ComplexMatrix operator*(ComplexMatrix::complex_type c, const ComplexMatrix& rhs)
{
    return rhs * c;
}

// Element wise division for matrix
ComplexMatrix ComplexMatrix::operator/(const ComplexMatrix& rhs) const
{
    assert(m_rows == rhs.getRowCount());
    assert(m_columns == rhs.getColumnCount());

    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] / rhs(rr, cc);
        }
    }
    return m;
}

ComplexMatrix& ComplexMatrix::operator/=(const ComplexMatrix& rhs)
{
    *this = *this / rhs;
    return *this;
}

ComplexMatrix ComplexMatrix::operator/(const complex_type& rhs) const
{
    ComplexMatrix m {m_rows, m_columns};

    for (unsigned rr=0; rr < m_rows; ++rr)
    {
        for (unsigned cc=0; cc < m_columns; ++cc)
        {
            m(rr, cc) = m_matrix[rr][cc] / rhs;
        }
    }
    return m;
}

ComplexMatrix operator/(ComplexMatrix::complex_type c, const ComplexMatrix& rhs)
{
    return rhs / c;
}

ComplexMatrix ComplexMatrix::subMatrix(unsigned rowStart, unsigned rowEnd, unsigned columnStart, unsigned columnEnd) const
{
    assert(rowStart < rowEnd);
    assert(columnStart < columnEnd);
    assert(rowStart < m_rows);
    assert(columnStart < m_columns);

    ComplexMatrix sub {rowEnd-rowStart+1, columnEnd-columnStart+1};

    for (unsigned rr=rowStart; rr <= rowEnd; ++rr)
    {
        for (unsigned cc=columnStart; cc <= columnEnd; ++cc)
        {
            sub(rr-rowStart, cc-columnStart) = m_matrix[rr][cc];
        }
    }

    return sub;
}

ComplexMatrix ComplexMatrix::transpose() const
{
  ComplexMatrix m {m_columns, m_rows};

  for (unsigned rr = 0; rr < m_columns; ++rr)
  {
    for (unsigned cc = 0; cc < m_rows; ++cc)
    {
      m(rr,cc) = m_matrix[cc][rr];
    }
  }

  return m;
}

ComplexMatrix ComplexMatrix::conjugateTranspose() const
{
  ComplexMatrix m {m_columns, m_rows};

  complex_type tmp;
  for (unsigned rr = 0; rr < m_columns; ++rr)
  {
    for (unsigned cc = 0; cc < m_rows; ++cc)
    {
        tmp = m_matrix[cc][rr];
        m(rr,cc) = complex_type{tmp.real(), -tmp.imag()};
    }
  }

  return m;
}


ComplexMatrix ComplexMatrix::minorMatrix(unsigned row, unsigned col) const
{
    assert(m_rows >= 2);
    assert(m_columns >= 2);
    assert(row < m_rows);
    assert(col < m_columns);

    ComplexMatrix newMatrix {m_rows - 1, m_columns - 1};

    for (unsigned rr=0; rr < m_rows-1; ++rr)
    {
        for (unsigned cc=0; cc < m_columns-1; ++cc)
        {
            newMatrix(rr, cc) = m_matrix[rr + (rr >= row)][cc + (cc >= col)];
        }
    }

  return newMatrix;
}

ComplexMatrix::complex_type ComplexMatrix::determinant() const
{
    assert(m_rows == m_columns);

    complex_type det {0.0, 0.0};

    if (m_rows == 1)
    {
      det = m_matrix[0][0];
    }
    else if (m_rows == 2)
    {
      // determinant of [a11,a12;a21,a22] is det = a11*a22-a21*a12
      det = m_matrix[0][0] * m_matrix[1][1] - m_matrix[1][0] * m_matrix[0][1];
    }
    else
    {
        double sign = 1.0;
        ComplexMatrix newMatrix (m_rows-1, m_columns-1);
        for (unsigned cc = 0; cc < m_columns; ++cc )
        {
            newMatrix = this->minorMatrix(0, cc);
//            det += pow(-1.0, cc) * m_matrix[0][cc] * newMatrix.determinant();
            det += sign * m_matrix[0][cc] * newMatrix.determinant();
            sign *= -1.0;
      }
    }

    return det;
}

void ComplexMatrix::swapValues(unsigned rowA, unsigned columnA, unsigned rowB, unsigned columnB)
{
    assert(rowA < m_rows);
    assert(rowB < m_rows);
    assert(columnA < m_columns);
    assert(columnB < m_columns);

    complex_type tmp {m_matrix[rowA][columnA]};
    m_matrix[rowA][columnA] = m_matrix[rowB][columnB];
    m_matrix[rowB][columnB] = tmp;
}

/*
 * returns the inverse of Matrix a
 */
ComplexMatrix ComplexMatrix::inverse()
{
    assert(m_rows == m_columns);

    ComplexMatrix res {m_rows, m_columns};

    if (m_rows == 1)
    {
        res(0, 0) = 1.0 / m_matrix[0][0];
    }
    else if (m_rows == 2)
    {
        complex_type d = this->determinant();
        assert( fabs(d) != 0 );
        res(0, 0) =  m_matrix[1][1];
        res(0, 1) = -m_matrix[0][1];
        res(1, 0) = -m_matrix[1][0];
        res(1, 1) =  m_matrix[0][0];
        res = res / d;
    }
    else
    {
        // 3 x 3 or larger
        // calculate inverse using gauss-jordan elimination
        //   http://mathworld.wolfram.com/MatrixInverse.html
        //   http://math.uww.edu/~mcfarlat/inverse.htm
        res = identity(m_rows);
        complex_type zero {0.0, 0.0};

        ComplexMatrix ai = *this; // make a copy of Matrix a

        for ( unsigned cc = 0; cc < m_columns; ++cc )
        {
            // element (cc, cc) should be non zero. if not, swap content
            // of lower rows
            unsigned rr;
            for ( rr = cc; (rr < m_rows) && (ai(rr, cc) == zero) ; ++rr )
            {
//                std::cout << cc << "  " << rr << " - " << ai(rr, cc) << std::endl;
            }
            assert(rr < m_rows); // If this fails, the entire row was zero
            if ( rr != cc )
            {
                // swap rows
                for (unsigned ss = 0; ss < m_columns ; ++ss )
                {
                    ai.swapValues(cc, ss, rr, ss);
                    res.swapValues(cc, ss, rr, ss);
                }
            }

            // eliminate non-zero values on the other rows at column cc
            for ( unsigned rr = 0; rr < m_rows ; ++rr )
            {
                if ( rr != cc )
                {
                    // eleminate value at cc and rr
                    if ( ai(rr, cc) != zero )
                    {
                        complex_type f = -ai(rr, cc) / ai(cc, cc);

                        // add (f * row cc) to row rr to eleminate the value
                        // at column cc
                        for ( unsigned ss = 0; ss < m_columns ; ++ss )
                        {
                            ai(rr, ss) += f * ai(cc, ss);
                            res(rr, ss) += f * res(cc, ss);
                        }
                    }
                }
                else
                {
                    // make value at (cc, cc) one,
                    // divide each value on row rr with the value at ai(cc, cc)
                    complex_type f = ai(cc, cc);
                    for ( unsigned ss = 0; ss < m_columns; ++ss )
                    {
                        ai(rr, ss) /= f;
                        res(rr, ss) /= f;
                    }
                }
            }
        }
    }

    return res;
}

unsigned ComplexMatrix::getRowCount() const
{
    return m_rows;
}

unsigned ComplexMatrix::getColumnCount() const
{
  return m_columns;
}

void ComplexMatrix::print() const
{
    for (unsigned rr = 0; rr < m_rows ; ++rr)
    {
        for (unsigned cc = 0; cc < m_columns ; ++cc)
        {
            std::cout << m_matrix[rr][cc] << " | ";
        }
        std::cout << std::endl;
    }
}
