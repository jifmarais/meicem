#include "LabelContainer.hpp"
#include <algorithm>
#include <assert.h>

LabelContainer::LabelContainer()
{
    //ctor
}

LabelContainer::~LabelContainer()
{
    //dtor
}

void LabelContainer::add(std::string label, SizeType index)
{
    // Add entry if it is not already in the map (no duplicates)

    std::map<std::string, std::vector<SizeType>>::const_iterator mapIterator;
    mapIterator = m_map.find(label);

    if ( mapIterator == m_map.cend() )
    {
        m_map[label].push_back(index);
    }
    else
    {
        const std::vector<SizeType>& indexList = mapIterator->second;
        if ( std::find(indexList.cbegin(), indexList.cend(), index) == indexList.cend() )
        {
            m_map[label].push_back(index);
        }
    }
}

std::vector<LabelContainer::SizeType> LabelContainer::find(std::string label) const
{
    std::vector<LabelContainer::SizeType> indexList;
    std::map<std::string, std::vector<SizeType>>::const_iterator mapIterator;
    mapIterator = m_map.find(label);

    if ( mapIterator != m_map.cend() )
    {
        indexList = mapIterator->second;
    }

    return indexList;
}

LabelContainer::SizeType LabelContainer::size() const
{
    return m_map.size();
}

