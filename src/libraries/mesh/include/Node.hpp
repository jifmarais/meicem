#pragma once

class Node
{
    public:

        Node();
        Node(double, double, double);
        virtual     ~Node();

        void        set(double, double, double);
        bool        operator==(const Node& rhs) const;
        bool        operator!=(const Node& rhs) const;
        double      x() const;
        double      y() const;
        double      z() const;

    protected:
    private:
        double m_x;
        double m_y;
        double m_z;
};

