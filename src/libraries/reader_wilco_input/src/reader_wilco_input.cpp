#include <fstream>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <boost/algorithm/string.hpp>
#include "Triangle.hpp"
#include "LabelContainer.hpp"
#include "reader_wilco_input.hpp"

WilcoInputReader::WilcoInputReader()
{
    m_TriangleContainer = nullptr;
    m_filename = "";
    m_frequency = 0.0;
}

WilcoInputReader::~WilcoInputReader()
{
    //dtor
}

void WilcoInputReader::setFile(std::string filename)
{
    m_filename = filename;
}

void WilcoInputReader::setTriangleContainer(TriangleContainer *triangles)
{
    m_TriangleContainer = triangles;
}

bool WilcoInputReader::hasTriangleContainer() const
{
    return m_TriangleContainer != nullptr;
}

void WilcoInputReader::unsetTriangleContainer()
{
   m_TriangleContainer = nullptr;
}

double WilcoInputReader::getFrequency() const
{
    return m_frequency;
}

bool WilcoInputReader::importModel() const
{
    assert(m_filename != "");
    std::ifstream file (m_filename);

    std::string line;
    std::vector<std::string> tokens;

    Node p;
    LabelContainer lContainer;


    // Using the triangle container's existing point container, if possible (all points are imported, used by triangles or not)
    NodeContainer* pContainer;
    if ( hasTriangleContainer() )
    {
        pContainer = &(m_TriangleContainer->getPointContainer());
    }
    else
    {
        NodeContainer internalNodeContainer;
        pContainer = &internalNodeContainer;
    }


    // Ensure that file exists
    assert(file);

    // Read the frequency
    getline(file,line);
    boost::algorithm::trim(line);
    m_frequency = std::stod(line);

    // Read the number of nodes
    getline(file,line);
    boost::algorithm::trim(line);
    int numNodes = std::stoi(line);

    // Read and process nodes
    int counter = -1;
    while (numNodes > 0)
    {
        numNodes -= 1;
        counter  += 1;
        tokens.clear();

        getline(file,line);
        boost::algorithm::trim(line);
        tokenizeLine(line, tokens);
        assert(tokens.size() >= 3);
//        printTokens(tokens);

        double x = std::stod(tokens.at(0));
        double y = std::stod(tokens.at(1));
        double z = std::stod(tokens.at(2));
        p.set(x, y, z);
        NodeContainer::SizeType index = pContainer->add(p);
        lContainer.add(std::to_string(counter), index);
    }

    // Read the number of triangles
    getline(file,line);
    boost::algorithm::trim(line);
    int numTriangles = std::stoi(line);

    // Read and process triangles
    if ( hasTriangleContainer() )
    {
        counter = 0;
        while (numTriangles > 0)
        {
            numTriangles -= 1;
            counter  += 1;
            tokens.clear();

            getline(file,line);
            boost::algorithm::trim(line);
            tokenizeLine(line, tokens);
            assert(tokens.size() >= 3);
//            printTokens(tokens);

            Triangle t;
            std::vector<NodeContainer::SizeType> p1IndexList = lContainer.find(tokens.at(0));
            std::vector<NodeContainer::SizeType> p2IndexList = lContainer.find(tokens.at(1));
            std::vector<NodeContainer::SizeType> p3IndexList = lContainer.find(tokens.at(2));
            assert(p1IndexList.size() == 1);
            assert(p2IndexList.size() == 1);
            assert(p3IndexList.size() == 1);
            t.set(pContainer->at(p1IndexList[0]),
                  pContainer->at(p2IndexList[0]),
                  pContainer->at(p3IndexList[0]));
            m_TriangleContainer->add(t);
        }
    }

    file.close();

    return false;
}

void WilcoInputReader::printTokens(const std::vector<std::string>& tokens) const
{
    int ii = 0;
    std::vector<std::string>::const_iterator token;
    for(token=tokens.begin() ; token < tokens.end(); ++token ) {
        std::string field = *token;
        std::cout << ii << "|" << field << "|" << std::endl;
        ii += 1;
    }
    std::cout << std::endl;
}

/* Tokenize the fields in the line based on the line formatting
  */
WilcoInputReader::returnResult WilcoInputReader::tokenizeLine(const std::string line,
                                                        std::vector<std::string> &tokens ) const
{
    std::string field;

    std::string::size_type lowerBoundary = 0;
    std::string::size_type upperBoundary = line.find(" ", lowerBoundary);

    while (upperBoundary != std::string::npos) {
        field = line.substr(lowerBoundary, upperBoundary - lowerBoundary);
        boost::algorithm::trim(field);
        tokens.push_back( field );
        lowerBoundary = upperBoundary + 1;
        upperBoundary = line.find(" ", lowerBoundary);
    }

    field = line.substr(lowerBoundary, upperBoundary - lowerBoundary);
    boost::algorithm::trim(field);
    tokens.push_back( field );

    return SUCCESS;
}
