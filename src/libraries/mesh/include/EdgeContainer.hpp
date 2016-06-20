#pragma once
#include <vector>
#include <limits>
#include "TriangleContainer.hpp"
#include "NodeContainer.hpp"
#include "Edge.hpp"

class EdgeContainer
{
    public:
        typedef std::vector<double>::size_type SizeType;
        static const SizeType invalidIndex = std::numeric_limits<SizeType>::max();

        EdgeContainer (TriangleContainer& container);
        virtual     		~EdgeContainer();

        TriangleContainer&      getTriangleContainer() const;
        void                    buildNonboundaryEdgeList();

        Edge                    at(SizeType index) const;
        SizeType    		size() const;

    protected:
    private:
//        SizeType    		add(Edge t);
//        bool                    isBoundaryEdge(SizeType index1) const;

        TriangleContainer&      m_triangleContainer;
        std::vector<NodeContainer::SizeType>	m_node1Index;
        std::vector<NodeContainer::SizeType>	m_node2Index;
        std::vector<std::vector<TriangleContainer::SizeType>> m_associatedTriangleIndeces;
};

