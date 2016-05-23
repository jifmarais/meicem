#include "Point3DContainer.hpp"
#include <algorithm>
#include <assert.h>
#include <cmath>
//#include <iostream>

Point3DContainer::Point3DContainer()
{
    m_size = 0;
    m_tolerance = 1e-6;
}

Point3DContainer::~Point3DContainer()
{
    //dtor
}

void Point3DContainer::setTolerance(double tolerance)
{
    assert(tolerance > 0);

    // We don't allow changing the tolerance once there are points in the list
    assert(m_size == 0);

    m_tolerance = tolerance;
}

double Point3DContainer::getTolerance() const
{
    return m_tolerance;
}

Point3DContainer::SizeType Point3DContainer::add(const Point3D &p)
{
    // Check if the point is already in the container
    SizeType index = find(p);
    if ( index == Point3DContainer::invalidIndex )
    {
        index = m_size;
        m_x.push_back(p.x());
        m_y.push_back(p.y());
        m_z.push_back(p.z());
        m_size += 1;
    }

    return index;
}

Point3DContainer::SizeType Point3DContainer::find(const Point3D &p) const
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
    return Point3DContainer::invalidIndex;    
}

Point3D Point3DContainer::at(SizeType index) const
{
    assert(index <= m_size);

    Point3D p;
    p.set(m_x[index], m_y[index], m_z[index]);
    return p;
}


Point3DContainer::SizeType Point3DContainer::size() const
{
    return m_size;
}

bool Point3DContainer::isEqual(double n1, double n2) const
{
    return std::abs(n1 - n2) <= m_tolerance;
}

