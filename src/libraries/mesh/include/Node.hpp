#pragma once

class Node
{
    public:

        Node();
        Node(double, double, double);
        virtual     ~Node();

        void        set(double x, double y, double z, double tolerance);
        void        set(double x, double y, double z);
        bool        operator==(const Node& rhs) const;
        bool        operator!=(const Node& rhs) const;
        Node&       operator=(const Node& rhs);
        Node        operator+(const Node& rhs) const;
        Node&       operator+=(const Node& rhs);
        Node        operator-(const Node& rhs) const;
        Node&       operator-=(const Node& rhs);
        Node        operator*(const double rhs) const;
        Node&       operator*=(const double rhs);
        Node        operator/(const double rhs) const;
        Node&       operator/=(const double rhs);
        double      magnitude() const;
        Node        norm() const;
        Node        cross(const Node& v) const;
        static Node cross(const Node& u, const Node& v);
        double 	    distance(const Node& p1) const;
        static double distance(const Node& p1, const Node& p2);
        double      x() const;
        double      y() const;
        double      z() const;

    protected:
    private:
        double m_x;
        double m_y;
        double m_z;
        double m_tolerance;

        bool   isEqual(double n1, double n2) const;
};

Node operator*(double d, const Node& rhs);
