#include "Triangle.hpp"

Triangle::Triangle()
{
    Point3D p {0.0, 0.0, 0.0};
    set(p , p, p);
}

Triangle::Triangle(Point3D n1,
                   Point3D n2,
                   Point3D n3)
{
   set(n1, n2, n3);
}

Triangle::~Triangle()
{
    //dtor
}

void Triangle::set(Point3D n1,
                   Point3D n2,
                   Point3D n3)
{
    m_nodes[0] = n1;
    m_nodes[1] = n2;
    m_nodes[2] = n3;
}

bool Triangle::operator==(const Triangle& rhs) const
{
    return ( m_nodes[0] == rhs.n1() &&
             m_nodes[1] == rhs.n2() &&
             m_nodes[2] == rhs.n3() );
}

bool Triangle::operator!=(const Triangle& rhs) const
{
    return not ( this->operator==(rhs) );
}

bool operator==(const Triangle& lhs, const Triangle& rhs)
{
    return lhs.operator==(rhs);
}

bool operator!=(const Triangle& lhs, const Triangle& rhs)
{
    return lhs.operator!=(rhs);
}

Point3D Triangle::n1() const
{
    return m_nodes[0];
}

Point3D Triangle::n2() const
{
    return m_nodes[1];
}

Point3D Triangle::n3() const
{
    return m_nodes[2];
}
