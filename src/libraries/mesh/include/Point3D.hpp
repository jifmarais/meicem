#pragma once

/**
 * @file Point3D.hpp
 * Write description of source file here for dOxygen.
 *
 * @brief Can use "brief" tag to explicitly generate comments for file documentation.
 *
 * @author Me 
 * @version 1.69 
 */

class Point3D
{
    public:

        Point3D();
        Point3D(double, double, double);
        virtual     ~Point3D();

        void        set(double, double, double);
        bool        operator==(const Point3D& rhs) const;
        bool        operator!=(const Point3D& rhs) const;
        double      x() const;
        double      y() const;
        double      z() const;

    protected:
    private:
        double m_x;
        double m_y;
        double m_z;
};

