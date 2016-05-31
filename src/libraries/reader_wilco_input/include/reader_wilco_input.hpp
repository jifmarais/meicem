#pragma once

#include <vector>
#include <string>
#include "TriangleContainer.hpp"

class WilcoInputReader
{
    public:

        WilcoInputReader();
        virtual		        ~WilcoInputReader();
        void 			setFile(std::string filename);
        void 			setTriangleContainer(TriangleContainer* triangles);
        bool			importModel() const;
        double 			getFrequency() const;

    protected:
    private:
        enum 		returnResult
                        {
                                SUCCESS,
                                FAIL,
                                NOT_IMPLEMENTED
                        };

        std::string		m_filename;
        TriangleContainer*	m_TriangleContainer;
        mutable double		m_frequency;

        bool			hasTriangleContainer() const;
        void			unsetTriangleContainer();

        std::vector<std::string> getNextLineOfTokens(std::ifstream &file) const;
        void 			printTokens(const std::vector<std::string>& tokens) const;
        returnResult 		tokenizeLine(std::string line,
                                             std::vector<std::string> &tokens ) const;

};

