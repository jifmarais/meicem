#include <iostream>
#include "EdgeContainer.hpp"
#include "NodeContainer.hpp"
#include "TriangleContainer.hpp"
#include <algorithm>
#include <assert.h>

EdgeContainer::EdgeContainer(const TriangleContainer &container)
   : m_triangleContainer(container)
{
    //ctor
}

EdgeContainer::~EdgeContainer()
{
    //dtor
}

const TriangleContainer& EdgeContainer::getTriangleContainer() const
{
    return m_triangleContainer;
}

Edge EdgeContainer::at(SizeType index) const
{
    assert(index < size());
    Edge e;
    auto pContainer = m_triangleContainer.getPointContainer();
    e.set(pContainer.at(m_node1Index.at(index)), pContainer.at(m_node2Index.at(index)));
    std::vector<TriangleContainer::SizeType> tIndeces = m_associatedTriangleIndeces.at(index);
    for (unsigned ii = 0; ii < tIndeces.size() ; ++ii)
    {
        e.associateTriangle(tIndeces.at(ii));
    }
    return e;
}

void EdgeContainer::buildNonboundaryEdgeList()
{
    NodeContainer pContainer = m_triangleContainer.getPointContainer();

    // Loop over all triangles and build a unique list of non-boundary edges.
    for (unsigned indexOuter = 0; indexOuter < m_triangleContainer.size() ; ++ indexOuter)
    {
        Triangle tOuter;
        tOuter = m_triangleContainer.at(indexOuter);

        for (unsigned indexInner = indexOuter + 1; indexInner < m_triangleContainer.size() ; ++ indexInner)
        {
            if ( m_triangleContainer.hasCommonNode(indexInner, indexOuter) )
            {
                 Triangle tInner;
                 tInner = m_triangleContainer.at(indexInner);

                 // Find common nodes
                 std::vector<NodeContainer::SizeType> edgePointList;
                 for (unsigned ii = 0; ii < 3; ++ii)
                 {
                     for (unsigned jj = 0; jj < 3; ++jj)
                     {
                         if (tOuter[ii] == tInner[jj])
                         {
                             edgePointList.push_back(pContainer.find(tOuter[ii]));
                         }
                     }
                 }

                 assert(edgePointList.size() < 3); // Three matching nodes = coincident triangle
                 if (edgePointList.size() == 2)
                 {
                     // Found a common edge
                     std::sort(edgePointList.begin(), edgePointList.end());
                     m_node1Index.push_back(edgePointList.at(0));
                     m_node2Index.push_back(edgePointList.at(1));
                     m_associatedTriangleIndeces.push_back({indexOuter, indexInner});
                 }
            }

        }

    }
}

NodeContainer::SizeType EdgeContainer::size() const
{
    return m_node1Index.size();
}

std::vector<EdgeContainer::SizeType> EdgeContainer::getEdgeIndecesOnTriangle(EdgeContainer::SizeType tIndex) const
{
    //JIF: This method does not have a test
    std::vector<SizeType> indexList;

    for (unsigned eIndex = 0; eIndex < m_associatedTriangleIndeces.size(); ++eIndex)
    {
        for (unsigned tAssociatedIndex = 0; tAssociatedIndex < m_associatedTriangleIndeces.at(eIndex).size(); ++tAssociatedIndex)
        {
            if (m_associatedTriangleIndeces.at(eIndex).at(tAssociatedIndex) == tIndex)
            {
                // edge is associated with the triangle
                indexList.push_back(eIndex);
            }
        }
    }

    return indexList;
}

//EdgeContainer::SizeType EdgeContainer::add(Edge e)
//{
//    return 1;
//}

//bool EdgeContainer::isBoundaryEdge(SizeType index1) const
//{
//    return false;
//}
