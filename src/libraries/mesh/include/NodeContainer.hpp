#pragma once
#include <vector>
#include <limits>
#include "Node.hpp"

class NodeContainer
{
    public:
        typedef std::vector<double>::size_type SizeType;
        static const SizeType invalidIndex = std::numeric_limits<SizeType>::max();

        NodeContainer();
        virtual     	~NodeContainer();

        void        	setTolerance(double);
        double      	getTolerance() const;
        SizeType    	add(Node);
        SizeType	find(const Node&) const;
        const Node&		at(SizeType) const;
        SizeType    	size() const;
        void        	clear();

    protected:
    private:
        std::vector<double>     m_x;
        std::vector<double>     m_y;
        std::vector<double>     m_z;
        std::vector<Node>		m_nodeList;
        SizeType                m_size;
        double                  m_tolerance;

        bool                    isEqual(double, double) const;
};

