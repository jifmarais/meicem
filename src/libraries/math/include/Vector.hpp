#pragma once

class Vector
{
    public:

        Vector();
        Vector(double, double, double);
        virtual         ~Vector();


        void            set(double, double, double);
        bool            operator==(const Vector& rhs) const;
        bool            operator!=(const Vector& rhs) const;
        Vector&         operator=(const Vector& rhs);
        Vector          operator+(const Vector& rhs) const;
        Vector&         operator+=(const Vector& rhs);
        Vector          operator-(const Vector& rhs) const;
        Vector&         operator-=(const Vector& rhs);
        Vector          operator*(const double rhs) const;
        Vector&         operator*=(const double rhs);
        Vector          operator/(const double rhs) const;
        Vector&         operator/=(const double rhs);
        double          magnitude() const;
        Vector          norm() const;
        double          x() const;
        double          y() const;
        double          z() const;

        static Vector   cross(const Vector& u, const Vector& v);

    protected:
    private:
        double m_x;
        double m_y;
        double m_z;
};

Vector operator*(double d, const Vector& rhs);
