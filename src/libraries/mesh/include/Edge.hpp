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

        void    set(const Node node1, const Node node2);
        void    associateTriangle(TriangleContainer::SizeType triangleIndex);
        void    setSortedAssociatedTriangles(const std::vector<TriangleContainer::SizeType> sortedUniqueTriangleList);
        const std::vector<TriangleContainer::SizeType> &getTriangles();
        double  length() const;
        Node	normal() const;
        void	correctOrientation(Triangle t);
        const Node &n1() const;
        const Node &n2() const;

    protected:
    private:
        Node    m_nodes [2];
        std::vector<TriangleContainer::SizeType> m_triangleIndeces;
};

