#include <assert.h>
#include <algorithm>
#include "Edge.hpp"

Edge::Edge()
{
    Node n1;
    Node n2;
    set(n1, n2);
}

Edge::Edge(Node node1, Node node2)
{
   set(node1, node2);
}

Edge::~Edge()
{
    //dtor
}

void Edge::set(Node node1,
               Node node2)
{
    m_nodes[0] = node1;
    m_nodes[1] = node2;
}

void Edge::associateTriangle(TriangleContainer::SizeType triangleIndex)
{
    assert(triangleIndex != TriangleContainer::invalidIndex);
    m_triangleIndeces.push_back(triangleIndex);
    std::sort(m_triangleIndeces.begin(), m_triangleIndeces.end());
    m_triangleIndeces.erase( std::unique( m_triangleIndeces.begin(), m_triangleIndeces.end() ), m_triangleIndeces.end() );
}

std::vector<TriangleContainer::SizeType> Edge::getTriangles()
{
    return m_triangleIndeces;
}

double Edge::length() const
{
    return Node::distance(m_nodes[0], m_nodes[1]);
}

Node Edge::n1() const
{
    return m_nodes[0];
}

Node Edge::n2() const
{
    return m_nodes[1];
}
