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

ComplexMatrix::~ComplexMatrix()
{
    //dtor
}

void ComplexMatrix::init(unsigned rows, unsigned columns, const complex_type initialValue)
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

ComplexMatrix ComplexMatrix::operator+(const complex_type rhs) const
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

ComplexMatrix ComplexMatrix::operator-(const complex_type rhs) const
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

ComplexMatrix ComplexMatrix::operator*(const complex_type rhs) const
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

ComplexMatrix ComplexMatrix::operator/(const complex_type rhs) const
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
        for (unsigned cc = 0; cc < m_columns; ++cc )
        {
            ComplexMatrix newMatrix = this->minorMatrix(0, cc);
//            det += pow(-1.0, cc) * m_matrix[0][cc] * newMatrix.determinant();
            det += sign * m_matrix[0][cc] * newMatrix.determinant();
            sign *= -1.0;
      }
    }

    return det;
}

///// NEED TO ADD TEST FOR ITEMS BELOW (EXCEPT FOR ROUND BRACKED OPERATOR)
//// swap two values
//void ComplexMatrix::swapValues(complex_type a, complex_type b)
//{
//  double tmp = a;
//  a = b;
//  b = tmp;
//}

///*
// * returns the inverse of Matrix a
// */
//Matrix Inv(const Matrix& a)
//{
//  Matrix res;
//  double d = 0;    // value of the determinant
//  int rows = a.GetRows();
//  int cols = a.GetRows();

//  d = Det(a);
//  if (rows == cols && d != 0)
//  {
//    // this is a square matrix
//    if (rows == 1)
//    {
//      // this is a 1 x 1 matrix
//      res = Matrix(rows, cols);
//      res(1, 1) = 1 / a.get(1, 1);
//    }
//    else if (rows == 2)
//    {
//      // this is a 2 x 2 matrix
//      res = Matrix(rows, cols);
//      res(1, 1) = a.get(2, 2);
//      res(1, 2) = -a.get(1, 2);
//      res(2, 1) = -a.get(2, 1);
//      res(2, 2) = a.get(1, 1);
//      res = (1/d) * res;
//    }
//    else
//    {
//      // this is a matrix of 3 x 3 or larger
//      // calculate inverse using gauss-jordan elimination
//      //   http://mathworld.wolfram.com/MatrixInverse.html
//      //   http://math.uww.edu/~mcfarlat/inverse.htm
//      res = Diag(rows);   // a diagonal matrix with ones at the diagonal
//      Matrix ai = a;    // make a copy of Matrix a

//      for (int c = 1; c <= cols; c++)
//      {
//        // element (c, c) should be non zero. if not, swap content
//        // of lower rows
//        int r;
//        for (r = c; r <= rows && ai(r, c) == 0; r++)
//        {
//        }
//        if (r != c)
//        {
//          // swap rows
//          for (int s = 1; s <= cols; s++)
//          {
//            Swap(ai(c, s), ai(r, s));
//            Swap(res(c, s), res(r, s));
//          }
//        }

//        // eliminate non-zero values on the other rows at column c
//        for (int r = 1; r <= rows; r++)
//        {
//          if(r != c)
//          {
//            // eleminate value at column c and row r
//            if (ai(r, c) != 0)
//            {
//              double f = - ai(r, c) / ai(c, c);

//              // add (f * row c) to row r to eleminate the value
//              // at column c
//              for (int s = 1; s <= cols; s++)
//              {
//                ai(r, s) += f * ai(c, s);
//                res(r, s) += f * res(c, s);
//              }
//            }
//          }
//          else
//          {
//            // make value at (c, c) one,
//            // divide each value on row r with the value at ai(c,c)
//            double f = ai(c, c);
//            for (int s = 1; s <= cols; s++)
//            {
//              ai(r, s) /= f;
//              res(r, s) /= f;
//            }
//          }
//        }
//      }
//    }
//  }
//  else
//  {
//    if (rows == cols)
//    {
//      throw Exception("Matrix must be square");
//    }
//    else
//    {
//      throw Exception("Determinant of matrix is zero");
//    }
//  }
//  return res;
//}

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

    for (unsigned rr = 0; rr < m_columns ; ++rr)
    {
        for (unsigned cc = 0; cc < m_rows ; ++cc)
        {
                std::cout << m_matrix[rr][cc] << " | ";
        }
        std::cout << std::endl;
    }

}
