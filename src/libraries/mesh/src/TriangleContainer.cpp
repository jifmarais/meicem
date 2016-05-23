#include "TriangleContainer.hpp"
#include <assert.h>

TriangleContainer::TriangleContainer(Point3DContainer& container)
   : m_pointContainer(container)
{
    m_tolerance = 1e-6;
}

TriangleContainer::~TriangleContainer()
{
    //dtor
}

Point3DContainer& TriangleContainer::getPointContainer() const
{
    return m_pointContainer;
}

TriangleContainer::SizeType TriangleContainer::add(const Triangle t)
{
    Point3DContainer::SizeType index;
    Point3DContainer::SizeType pIndex;
    index = m_node1.size();

    pIndex = m_pointContainer.add(t.n1());
    m_node1.push_back(pIndex);
    pIndex = m_pointContainer.add(t.n2());
    m_node2.push_back(pIndex);
    pIndex = m_pointContainer.add(t.n3());
    m_node3.push_back(pIndex);
    return index;
}

TriangleContainer::SizeType TriangleContainer::find(const TriangleContainer::SizeType i1,
                                                    const TriangleContainer::SizeType i2,
                                                    const TriangleContainer::SizeType i3) const
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

TriangleContainer::SizeType TriangleContainer::find(const Triangle t) const
{
    return find(m_pointContainer.find(t.n1()),
                m_pointContainer.find(t.n2()),
                m_pointContainer.find(t.n3()));
}

Triangle TriangleContainer::getAt(const SizeType index) const
{
    Point3D p1 = m_pointContainer.at(m_node1[index]);
    Point3D p2 = m_pointContainer.at(m_node2[index]);
    Point3D p3 = m_pointContainer.at(m_node3[index]);
    return Triangle {p1, p2, p3};
}

Point3DContainer::SizeType TriangleContainer::size() const
{
    return m_node1.size();
}

TriangleContainer::SizeType TriangleContainer::matchIndices(TriangleContainer::SizeType ii,
                                                            TriangleContainer::SizeType i1,
                                                            TriangleContainer::SizeType i2,
                                                            TriangleContainer::SizeType i3) const
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
