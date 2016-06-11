#include "Triangle.hpp"

Triangle::Triangle()
{
    Node p {0.0, 0.0, 0.0};
    set(p , p, p);
}

Triangle::Triangle(Node n1,
                   Node n2,
                   Node n3)
{
   set(n1, n2, n3);
}

Triangle::~Triangle()
{
    //dtor
}

void Triangle::set(Node n1,
                   Node n2,
                   Node n3)
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

double Triangle::area() const
{
    return 0.5*Node::cross(m_nodes[1]-m_nodes[0], m_nodes[2]-m_nodes[0]).magnitude();
}

Node Triangle::centre() const
{
    return ( m_nodes[0] + m_nodes[1] + m_nodes[2] ) / 3.0;
}

Node Triangle::n1() const
{
    return m_nodes[0];
}

Node Triangle::n2() const
{
    return m_nodes[1];
}

Node Triangle::n3() const
{
    return m_nodes[2];
}
