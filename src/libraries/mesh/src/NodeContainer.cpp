#include "NodeContainer.hpp"
#include <algorithm>
#include <assert.h>
#include <cmath>
//#include <iostream>

NodeContainer::NodeContainer()
{
    m_size = 0;
    m_tolerance = 1e-6;
}

NodeContainer::~NodeContainer()
{
    //dtor
}

void NodeContainer::setTolerance(double tolerance)
{
    assert(tolerance > 0);

    // We don't allow changing the tolerance once there are points in the list
    assert(m_size == 0);

    m_tolerance = tolerance;
}

double NodeContainer::getTolerance() const
{
    return m_tolerance;
}

NodeContainer::SizeType NodeContainer::add(const Node &p)
{
    // Check if the point is already in the container
    SizeType index = find(p);
    if ( index == NodeContainer::invalidIndex )
    {
        index = m_size;
        m_x.push_back(p.x());
        m_y.push_back(p.y());
        m_z.push_back(p.z());
        m_size += 1;
    }

    return index;
}

NodeContainer::SizeType NodeContainer::find(const Node &p) const
{
    // CRC:JIF: Would it be better to try to use std::find for this instead of a FOR loop?

    for ( SizeType ii = 0; ii < m_size ; ++ii )
    {
        if ( isEqual(p.x(), m_x[ii]) )
        {
            if ( isEqual(p.y(), m_y[ii]) )
            {
                if ( isEqual(p.z(), m_z[ii]) )
                {
                    return ii;
                }
            }
        }
    }
    return NodeContainer::invalidIndex;
}

Node NodeContainer::at(SizeType index) const
{
    assert(index < size());

    Node p;
    p.set(m_x[index], m_y[index], m_z[index], m_tolerance);
    return p;
}


NodeContainer::SizeType NodeContainer::size() const
{
    return m_size;
}

bool NodeContainer::isEqual(double n1, double n2) const
{
    return std::fabs(n1 - n2) <= m_tolerance;
}

void NodeContainer::clear()
{
    m_x.clear();
    m_y.clear();
    m_z.clear();
    m_size = 0;
}
