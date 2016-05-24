#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <assert.h>
#include <boost/algorithm/string.hpp>
#include "Triangle.hpp"
#include "LabelContainer.hpp"
#include "reader_nastran.hpp"

NastranReader::NastranReader()
{
    //ctor
}

NastranReader::~NastranReader()
{
    //dtor
}

void NastranReader::setFile(std::string filename)
{
    m_filename = filename;
}

void NastranReader::setTriangleContainer(TriangleContainer *triangles)
{
    m_TriangleContainer = triangles;
}

bool NastranReader::hasTriangleContainer() const
{
    return m_TriangleContainer != nullptr;
}

void NastranReader::unsetTriangleContainer()
{
   m_TriangleContainer = nullptr;
}

std::vector<std::string> NastranReader::getNextLineOfTokens(std::ifstream &file) const
{
    std::string line;
    std::vector<std::string> tokens;

    // Ensure that file exists
    assert(file);

    while ( getline(file,line) )
    {
        removeComments( line );
        lineFormat format;
        determineLineFormat(line, format);
        tokenizeLine(line, format, tokens);

        // Ensure that each line has 10 tokens
        if ( tokens.size() > 0 )
        {
            while ( tokens.size() < 10 )
            {
                tokens.push_back("");
            }
            break;
        }
    }

    return tokens;
}

bool NastranReader::importModel() const
{
    /*
      NASTRAN Statement (optional)
      File management statements (optional)
      Executive control statements (required)
      CEND
      Case Control Commands (required)
      BEGIN BULK
      Bulk Data Entries (required)
      ENDDATA
    */

    assert(m_filename != "");
    std::ifstream file (m_filename);

    std::string line;
    std::vector<std::string> tokens1;
    std::vector<std::string> tokens2;

    Node p;
    NodeContainer pContainer;
    LabelContainer lContainer;

    // Ensure that file exists
    assert(file);

    // Process file contents
    while (true)
    {
        tokens1 = getNextLineOfTokens(file);

        if ( tokens1.size() > 0 )
        {
            if ( tokens1.at(0) == "GRID" )
            {
//                printTokens(tokens1);
                double x = nasStringToDouble(tokens1.at(3));
                double y = nasStringToDouble(tokens1.at(4));
                double z = nasStringToDouble(tokens1.at(5));
                p.set(x, y, z);
                NodeContainer::SizeType index = pContainer.add(p);
                lContainer.add(tokens1.at(1), index);
            }
            else if ( tokens1.at(0) == "GRID*" )
            {
                tokens2 = getNextLineOfTokens(file);
                double x = nasStringToDouble(tokens1.at(3));
                double y = nasStringToDouble(tokens1.at(4));
                double z = nasStringToDouble(tokens2.at(1));
                p.set(x, y, z);
                NodeContainer::SizeType index = pContainer.add(p);
                lContainer.add(tokens1.at(1), index);
            }
        }
        else
        {
            break;
        }
    }

    // Reset input buffer to the start of the file
    file.clear();
    file.seekg(0, std::ios::beg);

    // Process file for elements
    while (true)
    {
        tokens1 = getNextLineOfTokens(file);

        if ( tokens1.size() > 0 )
        {
            if( tokens1.at(0) == "CTRIA3" )
            {
//                printTokens(tokens1);

                if ( hasTriangleContainer() )
                {
                    Triangle t;
                    std::vector<NodeContainer::SizeType> p1IndexList = lContainer.find(tokens1.at(3));
                    std::vector<NodeContainer::SizeType> p2IndexList = lContainer.find(tokens1.at(4));
                    std::vector<NodeContainer::SizeType> p3IndexList = lContainer.find(tokens1.at(5));
                    assert(p1IndexList.size() == 1);
                    assert(p2IndexList.size() == 1);
                    assert(p3IndexList.size() == 1);
                    t.set(pContainer.at(p1IndexList[0]),
                          pContainer.at(p2IndexList[0]),
                          pContainer.at(p3IndexList[0]));
                    m_TriangleContainer->add(t);
                }
            }
            else if( tokens1.at(0) == "CTRIA3*" )
            {
                // TODO
//                tokens2 = getNextLineOfTokens(file);

//                std::cout << tokens1.at(0) << std::endl;
//                std::cout << tokens2.at(0) << std::endl;
            }
            else
            {
//                std::cout << "Unhandled token: " << tokens1.at(0) << std::endl;
            }
        }
        else
        {
            break;
        }
    }

    file.close();

    return false;
}

void NastranReader::printTokens(const std::vector<std::string>& tokens) const
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

void NastranReader::removeComments(std::string &line) const
{
    std::string::size_type commentStartPossition = line.find("$");
    if (commentStartPossition != std::string::npos){
        line = line.substr(0,commentStartPossition);
    }
}

/*
     Theo following are valid decimal numbers in NASTRAN format that need to be supported
     7.0
     .70+1
     .7E1
     7.E+0
     0.7+1
     70.-1

     The numbers above can also be preceded with a + or a -.
   */
double NastranReader::nasStringToDouble(std::string token) const
{
    double real;

    std::string manipulatedToken = token;
    std::string significandString;
    std::string exponentString;

    double significand = 1;
    int exponent = 1;

    // If the token starts with + or -, remove it, but store its value
    if ( manipulatedToken.at(0) == '+' )
    {
        manipulatedToken = manipulatedToken.substr(1, manipulatedToken.size());
    }
    else if ( manipulatedToken.at(0) == '-' )
    {
        manipulatedToken = manipulatedToken.substr(1, manipulatedToken.size());
        significand = -1;
    }

    // If the token starts with ., add the missing zero
    if ( manipulatedToken.at(0) == '.' )
    {
        manipulatedToken.insert(manipulatedToken.cbegin(), '0');
    }

    // Now we can split the significand from the exponent at +, - or E
    std::string::size_type foundSplit = manipulatedToken.find_first_of("+-E");
    if ( (foundSplit != std::string::npos) )
    {
        significandString = manipulatedToken.substr(0, foundSplit);
        exponentString = manipulatedToken.substr(foundSplit, manipulatedToken.length() - foundSplit);

        // Remove any possible E from the exponent
        if ( exponentString.at(0) == 'E' )
        {
            exponentString = exponentString.substr(1, exponentString.size());
        }
    }
    else
    {
        significandString = manipulatedToken;
        exponentString = "0";
    }

    significand = significand * std::stod(significandString);
    exponent = exponent * std::stoi(exponentString);
    real = significand * pow( 10, exponent );

    return real;
}

/* Determine the format of the NASTRAN line: small, large, or free.
     Details from pyNastran - https://pynastran-locr.readthedocs.org/en/latest/manual/bdf_developer.html
     small: each line has 9 fields of 8 characters (72 chars max)
     large: 1x8, 4x16. Field 1 must have an asterisk following the character string
     free: each line entry must be separated by a comma
  */
NastranReader::returnResult NastranReader::determineLineFormat(std::string line, lineFormat & format ) const
{
    std::string::size_type foundAsterisk = line.find("*");
    std::string::size_type foundComma = line.find(",");

    if (foundComma == std::string::npos) {
        if(foundAsterisk == std::string::npos) {
            format = SMALL;
            return SUCCESS;
        } else {
            format = LARGE;
            return SUCCESS;
        }
    } else {
        if(foundAsterisk == std::string::npos) {
            format = SMALL_CSV;
            return SUCCESS;
        } else {
            format = LARGE_CSV;
            return SUCCESS;
        }
    }

    return FAIL;
}

/* Tokenize the fields in the line based on the line formatting
  */
NastranReader::returnResult NastranReader::tokenizeLine(const std::string line,
                                                        const lineFormat format,
                                                        std::vector<std::string> &tokens ) const
{
    std::string field;

    switch(format) {
    case SMALL: {
        // small: each line has maximum 10 fields of 8 characters
        const int fieldLength [10] {8, 8, 8, 8, 8, 8, 8, 8, 8, 8};

        std::string::size_type lowerBoundary = 0;
        for (unsigned int i=0; line.size() > lowerBoundary ; ++i) {
            field = line.substr( lowerBoundary,fieldLength[i] );
            boost::algorithm::trim(field);
            tokens.push_back( field );
            lowerBoundary += fieldLength[i];
        }
        break;
    }
    case LARGE: {
        // large: 1x8, 4x16, 1x8.
        const int fieldLength [6] {8, 16, 16, 16, 16, 8};

        std::string::size_type lowerBoundary = 0;
        for (unsigned int i=0; line.size() > lowerBoundary ; ++i) {
            field = line.substr( lowerBoundary,fieldLength[i] );
            boost::algorithm::trim(field);
            tokens.push_back( field );
            lowerBoundary += fieldLength[i];
        }
        break;
    }
    case SMALL_CSV:
    case LARGE_CSV: {
        std::string::size_type lowerBoundary = 0;
        std::string::size_type upperBoundary = line.find(",", lowerBoundary);

        while (upperBoundary != std::string::npos) {
            field = line.substr(lowerBoundary, upperBoundary - lowerBoundary);
            boost::algorithm::trim(field);
            tokens.push_back( field );
            lowerBoundary = upperBoundary + 1;
            upperBoundary = line.find(",", lowerBoundary);
        }

        field = line.substr(lowerBoundary, upperBoundary - lowerBoundary);
        boost::algorithm::trim(field);
        tokens.push_back( field );
        break;
    }
    default:
        return FAIL;
    }

    return SUCCESS;
}
