
//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

#include "reader_wilco_input.hpp"
#include "NodeContainer.hpp"
#include "TriangleContainer.hpp"

BOOST_AUTO_TEST_SUITE(WilcoInputFileReader_BasicTests)

std::string baseTestFilesDirectory = "../src/libraries/reader_wilco_input/test/test_files/";

BOOST_AUTO_TEST_CASE(noTriangleContainer)
{
    NodeContainer pContainer;
    WilcoInputReader reader;
    reader.setFile(baseTestFilesDirectory + "input.txt");
    reader.importModel();
    // Simply should not give an error
}

BOOST_AUTO_TEST_CASE(readFrequency)
{
    NodeContainer pContainer;
    WilcoInputReader reader;
    reader.setFile(baseTestFilesDirectory + "input.txt");
    reader.importModel();
    BOOST_CHECK_MESSAGE(reader.getFrequency() == 2400000000.0 , "The frequency value was not read correctly.");
}

BOOST_AUTO_TEST_CASE(readExampleFile)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    WilcoInputReader reader;
    reader.setFile(baseTestFilesDirectory + "input.txt");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 124, "There should be 124 triangles.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

