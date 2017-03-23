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

Node EdgeContainer::node1At(SizeType index) const
{
    assert(index < size());
    return m_edgeList.at(index).n1();
}

Node EdgeContainer::node2At(SizeType index) const
{
    assert(index < size());
    return m_edgeList.at(index).n2();
}

const std::vector<TriangleContainer::SizeType>& EdgeContainer::associatedTriaglesAt(SizeType index) const
{
    assert(index < size());
    return (m_edgeToTriangleIndecesMap.at(index));
}

const Edge& EdgeContainer::at(SizeType index) const
{
    assert(index < size());
    return m_edgeList.at(index);
}

void EdgeContainer::buildTriangleToEdgeMap()
{
    m_triangleToEdgeIndecesMap.clear();

    for (unsigned triangleIndex = 0; triangleIndex < m_triangleContainer.size(); ++ triangleIndex )
    {
        std::vector<EdgeContainer::SizeType> edgeList;
        for (unsigned edgeIndex = 0; edgeIndex < m_edgeToTriangleIndecesMap.size(); ++edgeIndex)
        {
            auto triangleList = m_edgeToTriangleIndecesMap.at(edgeIndex);
            for (unsigned assosiatedTriangleIndex = 0; assosiatedTriangleIndex < triangleList.size(); ++assosiatedTriangleIndex)
            {
                if (triangleList.at(assosiatedTriangleIndex) == triangleIndex)
                {
                    edgeList.push_back(edgeIndex);
                }
            }
        }
        m_triangleToEdgeIndecesMap.push_back(edgeList);
    }
}

void EdgeContainer::buildNonboundaryEdgeList()
{
    NodeContainer pContainer = m_triangleContainer.getPointContainer();
    m_edgeToTriangleIndecesMap.clear();
    m_edgeList.clear();

    // Loop over all triangles and build a unique list of non-boundary edges.
    for (unsigned indexOuter = 0; indexOuter < m_triangleContainer.size() ; ++ indexOuter)
    {
        const Triangle& tOuter = m_triangleContainer.at(indexOuter);

        for (unsigned indexInner = indexOuter + 1; indexInner < m_triangleContainer.size() ; ++ indexInner)
        {
            if ( m_triangleContainer.hasCommonNode(indexInner, indexOuter) )
            {
                 const Triangle& tInner = m_triangleContainer.at(indexInner);

                 // Find common nodes
                 std::vector<NodeContainer::SizeType> edgePointList;
                 for (unsigned ii = 0; ii < 3; ++ii)
                 {
                     for (unsigned jj = 0; jj < 3; ++jj)
                     {
                         if (tOuter.at(ii) == tInner.at(jj))
                         {
                             edgePointList.push_back(pContainer.find(tOuter.at(ii)));
                         }
                     }
                 }

                 assert(edgePointList.size() < 3); // Three matching nodes = coincident triangle
                 if (edgePointList.size() == 2)
                 {
                     // Found a common edge
                     Edge e;
                     std::sort(edgePointList.begin(), edgePointList.end());
                     e.set(pContainer.at(edgePointList.at(0)), pContainer.at(edgePointList.at(1)));
                     // Add items in sorted order (smallest to largest - basis function direction)
                     if (indexOuter < indexInner)
                     {
                         e.setSortedAssociatedTriangles({indexOuter, indexInner});
                         m_edgeToTriangleIndecesMap.push_back({indexOuter, indexInner});
                     }
                     else
                     {
                         e.setSortedAssociatedTriangles({indexInner, indexOuter});
                         m_edgeToTriangleIndecesMap.push_back({indexInner, indexOuter});
                     }
                     m_edgeList.push_back(e);
                 }
            }

        }

    }
    buildTriangleToEdgeMap();
}

NodeContainer::SizeType EdgeContainer::size() const
{
    return m_edgeList.size();
}

const std::vector<EdgeContainer::SizeType>& EdgeContainer::getEdgeIndecesOnTriangle(EdgeContainer::SizeType tIndex) const
{
    return m_triangleToEdgeIndecesMap.at(tIndex);
}
