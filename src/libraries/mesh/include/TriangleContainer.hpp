#pragma once
#include <vector>
#include <limits>
#include "NodeContainer.hpp"
#include "Triangle.hpp"

class TriangleContainer
{
    public:
        typedef std::vector<double>::size_type SizeType;
        static const SizeType invalidIndex = std::numeric_limits<SizeType>::max();

        TriangleContainer(NodeContainer& container);
        virtual     		~TriangleContainer();

        NodeContainer&		getPointContainer() const;
        SizeType    		add(Triangle t);
        SizeType    		find(Triangle t) const;
        SizeType                find(SizeType i1, SizeType i2, SizeType i3) const;
        const Triangle &at(SizeType index) const;
        SizeType    		size() const;
        bool                    hasCommonNode(SizeType index1, SizeType index2) const;

    protected:
    private:
        NodeContainer& 			m_pointContainer;
        std::vector<SizeType>	m_node1;
        std::vector<SizeType>	m_node2;
        std::vector<SizeType>	m_node3;
        std::vector<Triangle>   m_triangleList;
        double                  m_tolerance;
        SizeType                matchIndices(SizeType ii, SizeType i1, SizeType i2, SizeType i3) const;
};

