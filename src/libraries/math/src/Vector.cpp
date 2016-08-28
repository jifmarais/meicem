#include "Vector.hpp"
#include <math.h>

Vector::Vector()
{
   set(0.0, 0.0, 0.0); 
}

Vector::Vector(double x, double y, double z)
{
   set(x, y, z); 
}

Vector::~Vector()
{
    //dtor
}

void Vector::set(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

bool Vector::operator==(const Vector& rhs) const
{
    return ( m_x == rhs.x() &&
             m_y == rhs.y() &&
             m_z == rhs.z() );
}

bool Vector::operator!=(const Vector& rhs) const
{
    return not ( this->operator==(rhs) );
}

//bool operator==(const Vector& lhs, const Vector& rhs)
//{
//    return lhs.operator==(rhs);
//}

//bool operator!=(const Vector& lhs, const Vector& rhs)
//{
//    return lhs.operator!=(rhs);
//}

Vector& Vector::operator=(const Vector& rhs)
{
    m_x = rhs.x();
    m_y = rhs.y();
    m_z = rhs.z();
    return *this;
}

Vector Vector::operator+(const Vector& rhs) const
{
    Vector v;
    v.set( m_x + rhs.x(), m_y + rhs.y(), m_z + rhs.z() );
    return v;
}

Vector& Vector::operator+=(const Vector& rhs)
{
    *this = *this + rhs;
    return *this;
}

Vector Vector::operator-(const Vector& rhs) const
{
    Vector v;
    v.set( m_x - rhs.x(), m_y - rhs.y(), m_z - rhs.z() );
    return v;
}

Vector& Vector::operator-=(const Vector& rhs)
{
    *this = *this - rhs;
    return *this;
}

Vector Vector::operator*(const double rhs) const
{
    Vector v;
    v.set( m_x * rhs, m_y * rhs, m_z * rhs );
    return v;
}

//Vector operator*(const Vector rhs, double d)
//{
//    Vector v;
//    v = rhs * d;
//    return v;
//}

Vector operator*(double d, const Vector& rhs)
{
    Vector v;
    v = rhs * d;
    return v;
}

Vector& Vector::operator*=(const double rhs)
{
    *this = *this * rhs;
    return *this;
}

Vector Vector::operator/(const double rhs) const
{
    Vector v;
    v.set( m_x / rhs, m_y / rhs, m_z / rhs );
    return v;
}

Vector& Vector::operator/=(const double rhs)
{
    *this = *this / rhs;
    return *this;
}


double Vector::magnitude() const
{
    return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
}

Vector Vector::norm() const
{
    Vector r;
    double vMag = this->magnitude();

    if (vMag > 0 )
    {
        r = *this / vMag;
    }
    else
    {
        r.set(0.0, 0.0, 0.0);
    }

    return r;
}

double Vector::dot(const Vector& u, const Vector& v) //JIF: Need to add test
{
    double n;
    n  = u.x()*v.x() + u.y()*v.y() + u.z()*v.z();
    return n;
}

Vector Vector::cross(const Vector& u, const Vector& v)
{
    Vector n;
    n.set(u.y()*v.z() - u.z()*v.y() ,
          u.z()*v.x() - u.x()*v.z() ,
          u.x()*v.y() - u.y()*v.x() );
    return n;
}

double Vector::x() const
{
    return m_x;
}

double Vector::y() const
{
    return m_y;
}

double Vector::z() const
{
    return m_z;
}

