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

Node Edge::normal() const
{
   //JIF : Has no test
    Node n;
    n = m_nodes[1] - m_nodes[0];
    n = n / n.magnitude();
    return n;
}

void Edge::correctOrientation(Triangle t)
{
    //JIF: Has no test

    assert(t.n1() == m_nodes[0] || t.n2() == m_nodes[0] || t.n3() == m_nodes[0]);
    assert(t.n1() == m_nodes[1] || t.n2() == m_nodes[1] || t.n3() == m_nodes[1]);

    for (NodeContainer::SizeType nIndex = 0; nIndex < 3; ++nIndex )
    {
        if (t[nIndex] == m_nodes[0])
        {
            if (t[(nIndex+1)%3] != m_nodes[1])
            {
                Node tmp = m_nodes[0];
                m_nodes[0] = m_nodes[1];
                m_nodes[1] = tmp;
            }
            break;
        }
    }
}

Node Edge::n1() const
{
    return m_nodes[0];
}

Node Edge::n2() const
{
    return m_nodes[1];
}
