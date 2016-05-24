#pragma once
#include "Node.hpp"

class Triangle
{
    public:

        Triangle();
        Triangle(Node, Node, Node);
        virtual  ~Triangle();

        void    set(Node, Node, Node);
        bool    operator==(const Triangle& rhs) const;
        bool    operator!=(const Triangle& rhs) const;
        Node	n1() const;
        Node	n2() const;
        Node    n3() const;

    protected:
    private:
        Node    m_nodes [3];
};

