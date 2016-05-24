#include "Node.hpp"

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

void Node::set(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

bool Node::operator==(const Node& rhs) const
{
    return ( m_x == rhs.x() &&
             m_y == rhs.y() &&
             m_z == rhs.z() );
}

bool Node::operator!=(const Node& rhs) const
{
    return not ( this->operator==(rhs) );
}

bool operator==(const Node& lhs, const Node& rhs)
{
    return lhs.operator==(rhs);
}

bool operator!=(const Node& lhs, const Node& rhs)
{
    return lhs.operator!=(rhs);
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

