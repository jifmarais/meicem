#pragma once
#include <string>
#include <vector>
#include <map>
#include <limits>       // std::numeric_limits

class LabelContainer
{
    public:
        typedef std::vector<double>::size_type SizeType;
        static const SizeType invalidIndex = std::numeric_limits<SizeType>::max();

        LabelContainer();
        virtual     		~LabelContainer();

        void			 add(const std::string label, SizeType index);
        std::vector<SizeType>	find(const std::string label);
        SizeType    		size() const;

    protected:
    private:
        std::map<std::string, std::vector<SizeType>>	m_map;
};

