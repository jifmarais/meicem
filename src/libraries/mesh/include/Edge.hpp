#pragma once
#include "TriangleContainer.hpp"
#include "Node.hpp"

class Edge
{
    public:

        Edge();
        Edge(Node node1,
             Node node2);
        virtual  ~Edge();

        void    set(Node node1,
                    Node node2);
        void    associateTriangle(TriangleContainer::SizeType triangleIndex);
        std::vector<TriangleContainer::SizeType>
                getTriangles();
        double  length() const;
        Node    n1() const;
        Node    n2() const;

    protected:
    private:
        Node    m_nodes [2];
        std::vector<TriangleContainer::SizeType> m_triangleIndeces;
};

