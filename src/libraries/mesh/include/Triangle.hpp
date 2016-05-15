#pragma once
#include "Point3D.hpp"

/**
 * @file Triangle.hpp
 * Write description of source file here for dOxygen.
 *
 * @brief Can use "brief" tag to explicitly generate comments for file documentation.
 *
 * @author Me 
 * @version 1.69 
 */

class Triangle
{
    public:

        Triangle();
        Triangle(Point3D,
                 Point3D,
                 Point3D);
        virtual     ~Triangle();

        void        set(Point3D,
                       Point3D,
                       Point3D);
        bool        operator==(const Triangle& rhs) const;
        bool        operator!=(const Triangle& rhs) const;
        Point3D     n1() const;
        Point3D     n2() const;
        Point3D     n3() const;

    protected:
    private:
        Point3D     m_nodes [3];
};

