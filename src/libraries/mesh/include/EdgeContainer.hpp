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

        EdgeContainer (const TriangleContainer& container);
        virtual     		~EdgeContainer();

        const TriangleContainer&      getTriangleContainer() const;
        void                    buildNonboundaryEdgeList();
        Edge                    at(SizeType index) const;
        SizeType    		size() const;
        const std::vector<SizeType> &getEdgeIndecesOnTriangle(SizeType tIndex) const;

        Node node1At(SizeType index) const;
        Node node2At(SizeType index) const;
        const std::vector<TriangleContainer::SizeType> &associatedTriaglesAt(SizeType index) const;
protected:
    private:
        void buildTriangleToEdgeMap();

        const TriangleContainer&                m_triangleContainer;
        std::vector<NodeContainer::SizeType>	m_node1Index;
        std::vector<NodeContainer::SizeType>	m_node2Index;
        std::vector<std::vector<TriangleContainer::SizeType>> m_edgeToTriangleIndecesMap;
        std::vector<std::vector<EdgeContainer::SizeType>> 	  m_triangleToEdgeIndecesMap;
};

