#pragma once

#include <vector>
#include <string>
#include "TriangleContainer.hpp"

class NastranReader
{
    public:

        NastranReader();
        virtual     		~NastranReader();
        void 			setFile(std::string filename);
        void 			setTriangleContainer(TriangleContainer* triangles);
        bool			importModel() const;

    protected:
    private:
        enum 		returnResult
                        {
                                SUCCESS,
                                FAIL,
                                NOT_IMPLEMENTED
                        };
        enum 		lineFormat
                        {
                                SMALL,
                                LARGE,
                                SMALL_CSV,
                                LARGE_CSV,
                                UNKNOWN
                        };

        std::string		m_filename;
        TriangleContainer*	m_TriangleContainer;

        bool			hasTriangleContainer() const;
        void			unsetTriangleContainer();

        void 			printTokens(const std::vector<std::string>& tokens) const;
        std::vector<std::string> getNextLineOfTokens(std::ifstream &file) const;
        void 			removeComments(std::string &line) const;
        double 			nasStringToDouble(std::string token) const;
        returnResult 		determineLineFormat(std::string line, lineFormat &format ) const;
        returnResult 		tokenizeLine(std::string line, lineFormat format,
                                             std::vector<std::string> &tokens ) const;

};

