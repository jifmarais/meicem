#include "LabelContainer.hpp"
#include <assert.h>

LabelContainer::LabelContainer()
{
    //ctor
}

LabelContainer::~LabelContainer()
{
    //dtor
}

void LabelContainer::add(const std::string label, SizeType index)
{
    m_map[label].push_back(index);
}

std::vector<LabelContainer::SizeType> LabelContainer::find(const std::string label)
{
    std::vector<LabelContainer::SizeType> indexList;
    indexList = m_map[label];
    return indexList;
}

LabelContainer::SizeType LabelContainer::size() const
{
    return m_map.size();
}

