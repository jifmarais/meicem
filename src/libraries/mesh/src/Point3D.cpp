#include "Point3D.hpp"

Point3D::Point3D()
{
   set(0.0, 0.0, 0.0); 
}

Point3D::Point3D(double x, double y, double z)
{
   set(x, y, z); 
}

Point3D::~Point3D()
{
    //dtor
}

void Point3D::set(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

bool Point3D::operator==(const Point3D& rhs) const
{
    return ( m_x == rhs.x() &&
             m_y == rhs.y() &&
             m_z == rhs.z() );
}

bool Point3D::operator!=(const Point3D& rhs) const
{
    return not ( this->operator==(rhs) );
}

bool operator==(const Point3D& lhs, const Point3D& rhs)
{
    return lhs.operator==(rhs);
}

bool operator!=(const Point3D& lhs, const Point3D& rhs)
{
    return lhs.operator!=(rhs);
}

double Point3D::x() const
{
    return m_x;
}

double Point3D::y() const
{
    return m_y;
}

double Point3D::z() const
{
    return m_z;
}

