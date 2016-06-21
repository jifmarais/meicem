#include "Node.hpp"
#include <math.h>
#include <cmath>

Node::Node()
{
   set(0.0, 0.0, 0.0); 
}

Node::Node(double x, double y, double z)
{
   set(x, y, z); 
}

Node::~Node()
{
    //dtor
}

void Node::set(double x, double y, double z, double tolerance)
{
    m_x = x;
    m_y = y;
    m_z = z;
    m_tolerance = tolerance;
}

void Node::set(double x, double y, double z)
{
   set(x, y, z, 1e-6);
}

bool Node::operator==(const Node& rhs) const
{
    return ( isEqual(m_x, rhs.x()) &&
             isEqual(m_y, rhs.y()) &&
             isEqual(m_z, rhs.z()) );
}

bool Node::operator!=(const Node& rhs) const
{
    return not ( this->operator==(rhs) );
}

Node& Node::operator=(const Node& rhs)
{
    m_x = rhs.x();
    m_y = rhs.y();
    m_z = rhs.z();
    return *this;
}

bool operator==(const Node& lhs, const Node& rhs)
{
    return lhs.operator==(rhs);
}

bool operator!=(const Node& lhs, const Node& rhs)
{
    return lhs.operator!=(rhs);
}

Node Node::operator+(const Node& rhs) const
{
    Node v;
    v.set( m_x + rhs.x(), m_y + rhs.y(), m_z + rhs.z() );
    return v;
}

Node& Node::operator+=(const Node& rhs)
{
    *this = *this + rhs;
    return *this;
}

Node Node::operator-(const Node& rhs) const
{
    Node v;
    v.set( m_x - rhs.x(), m_y - rhs.y(), m_z - rhs.z() );
    return v;
}

Node& Node::operator-=(const Node& rhs)
{
    *this = *this - rhs;
    return *this;
}

Node Node::operator*(const double rhs) const
{
    Node v;
    v.set( m_x * rhs, m_y * rhs, m_z * rhs );
    return v;
}

//Node operator*(const Node rhs, double d)
//{
//    Node v;
//    v = rhs * d;
//    return v;
//}

Node operator*(double d, const Node& rhs)
{
    Node v;
    v = rhs * d;
    return v;
}

Node& Node::operator*=(const double rhs)
{
    *this = *this * rhs;
    return *this;
}

Node Node::operator/(const double rhs) const
{
    Node v;
    v.set( m_x / rhs, m_y / rhs, m_z / rhs );
    return v;
}

Node& Node::operator/=(const double rhs)
{
    *this = *this / rhs;
    return *this;
}


double Node::magnitude() const
{
    return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
}

Node Node::norm() const
{
    Node r;
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

Node Node::cross(const Node& u, const Node& v)
{
    Node n;
    n.set(u.y()*v.z() - u.z()*v.y() ,
          u.z()*v.x() - u.x()*v.z() ,
          u.x()*v.y() - u.y()*v.x() );
    return n;
}

Node Node::cross(const Node& v) const
{
    return cross(*this, v);
}

double Node::distance(const Node& p1, const Node& p2)
{
    return (p2 - p1).magnitude();
}

double Node::distance(const Node& p1) const
{
    return distance(*this, p1);
}


double Node::x() const
{
    return m_x;
}

double Node::y() const
{
    return m_y;
}

double Node::z() const
{
    return m_z;
}

bool Node::isEqual(double n1, double n2) const
{
    return std::fabs(n1 - n2) <= m_tolerance;
}
