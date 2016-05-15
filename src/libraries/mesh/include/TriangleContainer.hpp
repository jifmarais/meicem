#pragma once
#include <vector>
#include <limits>       // std::numeric_limits
#include "Point3DContainer.hpp"
#include "Triangle.hpp"

class TriangleContainer
{
    public:
        typedef std::vector<double>::size_type SizeType;
        static const SizeType invalidIndex = std::numeric_limits<SizeType>::max();

        TriangleContainer(Point3DContainer& container);
        virtual     ~TriangleContainer();

        Point3DContainer&       getPointContainer() const;
        SizeType    		add(const Triangle t);
        SizeType    		find(const Triangle t) const;
        SizeType                find(const TriangleContainer::SizeType i1,
                                     const TriangleContainer::SizeType i2,
                                     const TriangleContainer::SizeType i3) const;
        Triangle     		getAt(const SizeType index) const;
        SizeType    		size() const;

    protected:
    private:
        Point3DContainer& 	m_pointContainer;
        std::vector<SizeType>	m_node1;
        std::vector<SizeType>	m_node2;
        std::vector<SizeType>	m_node3;
        double                  m_tolerance;
        SizeType                matchIndices(TriangleContainer::SizeType ii,
                                             TriangleContainer::SizeType i1,
                                             TriangleContainer::SizeType i2,
                                             TriangleContainer::SizeType i3) const;
};

