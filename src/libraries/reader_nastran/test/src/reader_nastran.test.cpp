
//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

#include "reader_nastran.hpp"
#include "NodeContainer.hpp"
#include "TriangleContainer.hpp"

BOOST_AUTO_TEST_SUITE(NastranReader_BasicTests)

std::string baseTestFilesDirectory = "../src/libraries/reader_nastran/test/test_files/";

BOOST_AUTO_TEST_CASE(testNoContainersShouldBeAllowed)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "simple_triangles_long.nas");
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 0, "There should be 0 triangles since the triangle container was not set.");
}

BOOST_AUTO_TEST_CASE(simple_triangles_long)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "simple_triangles_long.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 36, "There should be 36 triangles in this import.");
}

BOOST_AUTO_TEST_CASE(ISat_Dploy_Sm)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "ISat_Dploy_Sm.dat");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 32, "There should be 32 triangles in this import.");
}

BOOST_AUTO_TEST_CASE(BWB_saero)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "BWB_saero.bdf");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 136, "There should be 136 triangles in this import.");
}

BOOST_AUTO_TEST_CASE(solid_shell_bar_xyz_bdf)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "solid_shell_bar_xyz.bdf");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 8, "There should be 8 triangles in this import.");
}

BOOST_AUTO_TEST_CASE(bend_A1_105)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "bend_A1_105.bdf");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 6, "There should be 6 triangles in this import.");
}

BOOST_AUTO_TEST_CASE(ctria3)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    reader.setFile(baseTestFilesDirectory + "ctria3.bdf");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    BOOST_CHECK_MESSAGE(tContainer.size() == 1, "There should be 1 triangle in this import.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

