#include <assert.h>
#include "Node.hpp"
#include <math.h>
#include <iostream>
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
    m_point(0) = x;
    m_point(1) = y;
    m_point(2) = z;
    m_tolerance = tolerance;
}

void Node::set(double x, double y, double z)
{
   set(x, y, z, 1e-6);
}

bool Node::operator==(const Node& rhs) const
{
    return ( isEqual(m_point(0), rhs.x()) &&
             isEqual(m_point(1), rhs.y()) &&
             isEqual(m_point(2), rhs.z()) );
}

bool Node::operator!=(const Node& rhs) const
{
    return not ( this->operator==(rhs) );
}

Node& Node::operator=(const Node& rhs)
{
    m_point(0) = rhs.x();
    m_point(1) = rhs.y();
    m_point(2) = rhs.z();
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
    v.set( m_point(0) + rhs.x(), m_point(1) + rhs.y(), m_point(2) + rhs.z() );
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
    v.set( m_point(0) - rhs.x(), m_point(1) - rhs.y(), m_point(2) - rhs.z() );
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
    v.set( m_point(0) * rhs, m_point(1) * rhs, m_point(2) * rhs );
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
    v.set( m_point(0) / rhs, m_point(1) / rhs, m_point(2) / rhs );
    return v;
}

Node& Node::operator/=(const double rhs)
{
    *this = *this / rhs;
    return *this;
}


double Node::magnitude() const
{
    return sqrt(m_point(0)*m_point(0) + m_point(1)*m_point(1) + m_point(2)*m_point(2));
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
    n.setVec(arma::cross(u.getVec(), v.getVec()));
//    n.set(u.y()*v.z() - u.z()*v.y() ,
//          u.z()*v.x() - u.x()*v.z() ,
//          u.x()*v.y() - u.y()*v.x() );
    return n;
}

double Node::dot(const Node &u, const Node &v)
{
    double n;
//    n  = u.x()*v.x() + u.y()*v.y() + u.z()*v.z();
    n = arma::dot(u.getVec(), v.getVec());
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

Node Node::transform(const arma::mat& transformMatrix) const
{
    assert(transformMatrix.n_rows == 3);
    assert(transformMatrix.n_cols == 3);

    Node newNode;
    newNode.setVec(transformMatrix*m_point);
//    newNode.set(transformMatrix(0,0)*m_point(0) + transformMatrix(0,1)*m_point(1) + transformMatrix(0,2)*m_point(2),
//                transformMatrix(1,0)*m_point(0) + transformMatrix(1,1)*m_point(1) + transformMatrix(1,2)*m_point(2),
//                transformMatrix(2,0)*m_point(0) + transformMatrix(2,1)*m_point(1) + transformMatrix(2,2)*m_point(2));
    return newNode;
}

void Node::setVec(vec p)
{
    m_point = p;
}

const Node::vec& Node::getVec() const
{
    return m_point;
}

double Node::x() const
{
    return m_point(0);
}

double Node::y() const
{
    return m_point(1);
}

double Node::z() const
{
    return m_point(2);
}

bool Node::isEqual(double n1, double n2) const
{
    return std::fabs(n1 - n2) <= m_tolerance;
}

void Node::print() const
{
    std::cout << "(" << m_point(0) << ", " << m_point(1) << ", " << m_point(2) << ")" << std::endl;
}
