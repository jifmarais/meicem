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
        bool			importModel();

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

        bool			hasTriangleContainer();
        void			unsetTriangleContainer();

        void 			printTokens(std::vector<std::string>& tokens);
        std::vector<std::string> getNextLineOfTokens(std::ifstream &file);
        void 			removeComments(std::string &line);
        double 			nasStringToDouble(std::string token);
        returnResult 		determineLineFormat(const std::string line,
                                        lineFormat &format );
        returnResult 		tokenizeLine(const std::string line,
                                        const lineFormat format,
                                        std::vector<std::string> &tokens );

};

