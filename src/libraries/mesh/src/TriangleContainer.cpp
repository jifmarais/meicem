#include "TriangleContainer.hpp"
#include <assert.h>

TriangleContainer::TriangleContainer(NodeContainer& container)
   : m_pointContainer(container)
{
    m_tolerance = 1e-6;
}

TriangleContainer::~TriangleContainer()
{
    //dtor
}

NodeContainer& TriangleContainer::getPointContainer() const
{
    return m_pointContainer;
}

TriangleContainer::SizeType TriangleContainer::add(Triangle t)
{
    NodeContainer::SizeType index;
    NodeContainer::SizeType pIndex;
    index = m_node1.size();

    pIndex = m_pointContainer.add(t.n1());
    m_node1.push_back(pIndex);
    pIndex = m_pointContainer.add(t.n2());
    m_node2.push_back(pIndex);
    pIndex = m_pointContainer.add(t.n3());
    m_node3.push_back(pIndex);

    m_triangleList.push_back(t);
    return index;
}

TriangleContainer::SizeType TriangleContainer::find(SizeType i1, SizeType i2, SizeType i3) const
{
    SizeType index = invalidIndex;
    for (SizeType ii = 0; ii < size() ; ++ii)
    {
        index = matchIndices(ii, i1, i2, i3);
        if ( index == invalidIndex )
        {
            index = matchIndices(ii, i3, i1, i2);
            if ( index == invalidIndex )
            {
                index = matchIndices(ii, i2, i3, i1);
                if ( index != invalidIndex )
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }
        else
        {
            break;
        }
    }
    return index;
}

TriangleContainer::SizeType TriangleContainer::find(Triangle t) const
{
    return find(m_pointContainer.find(t.n1()),
                m_pointContainer.find(t.n2()),
                m_pointContainer.find(t.n3()));
}

const Triangle& TriangleContainer::at(SizeType index) const
{
    assert(index < size());
    return m_triangleList.at(index);
}

NodeContainer::SizeType TriangleContainer::size() const
{
    return m_node1.size();
}

bool TriangleContainer::hasCommonNode(SizeType index1, SizeType index2) const
{
    Triangle t1 = at(index1);
    Triangle t2 = at(index2);

    for (unsigned ii = 0 ; ii < 3 ; ++ii )
    {
        for (unsigned jj = 0 ; jj < 3 ; ++jj )
        {
            if ( t1[ii] == t2[jj] )
            {
                return true;
            }
        }
    }
    return false;
}

TriangleContainer::SizeType TriangleContainer::matchIndices(SizeType ii,
                                                            SizeType i1,
                                                            SizeType i2,
                                                            SizeType i3) const
{
    SizeType index = invalidIndex;
    if ( i1 == m_node1[ii] )
    {
        if ( i2 == m_node2[ii] )
        {
            if ( i3 == m_node3[ii] )
            {
                index = ii;
            }
        }
        //else
        //{
            if ( i3 == m_node2[ii] )
            {
                if ( i2 == m_node3[ii] )
                {
                    index = ii;
                }
            }
        //}
    }
    return index;
}
