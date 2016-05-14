#pragma once

#include <vector>
#include <limits>       // std::numeric_limits
#include "Point3D.hpp"

/**
 * @file Point3DContainer.hpp
 * Write description of source file here for dOxygen.
 *
 * @brief Can use "brief" tag to explicitly generate comments for file documentation.
 *
 * @author Me 
 * @version 1.69 
 */

class Point3DContainer
{
    public:
        typedef std::vector<double>::size_type SizeType;
        static const SizeType invalidIndex = std::numeric_limits<SizeType>::max();

        Point3DContainer();
        virtual     ~Point3DContainer();

        void        setTolerance(double);
        double      getTolerance() const;
        SizeType    addPoint(const Point3D&);
        SizeType    findPoint(const Point3D&) const;
        Point3D     getPointAt(const SizeType) const;
        SizeType    size() const;

    protected:
    private:
        std::vector<double>     m_x;
        std::vector<double>     m_y;
        std::vector<double>     m_z;
        SizeType                m_size;
        double                  m_tolerance;

        bool                    isEqual(double, double) const;
};

