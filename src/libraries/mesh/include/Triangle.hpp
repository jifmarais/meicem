#pragma once
#include "Node.hpp"
#include "ComplexMatrix.hpp"

class Triangle
{
    public:

        Triangle();
        Triangle(Node, Node, Node);
        virtual  ~Triangle();

        void    set(Node, Node, Node);
        bool    operator==(const Triangle& rhs) const;
        bool    operator!=(const Triangle& rhs) const;
        Node    operator[](unsigned index) const;
        double  area() const;
        Node	centre() const;
        void    setOppositeEdge(const Node p1, const Node p2); // Needs tests
        Node	normal() const;
        Node	toSimplex(const Node& p) const;
        Node	fromSimplex(const Node& p) const;
        Node	n1() const;
        Node	n2() const;
        Node    n3() const;

        Triangle transform(const ComplexMatrix& transformMatrix) const; // Needs tests
protected:
    private:
        Node    m_nodes [3];
};

