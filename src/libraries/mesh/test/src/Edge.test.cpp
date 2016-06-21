#include <iostream>
#include "Edge.hpp"
#include "Node.hpp"
#include "TriangleContainer.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Edge_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    Node n1 {0.0, 0.0, 0.0};
    Node n2 {1.0, 1.0, 1.0};
    Edge e1 {n1, n2};
    BOOST_CHECK_MESSAGE(e1.n1() == n1, "Node not set correctly.");
    BOOST_CHECK_MESSAGE(e1.n2() == n2, "Node not set correctly.");

    Node n3 {2.0, 2.0, 2.0};
    e1.set(n3, n1);
    BOOST_CHECK_MESSAGE(e1.n1() == n3, "Node not set correctly.");
    BOOST_CHECK_MESSAGE(e1.n2() == n1, "Node not set correctly.");

    std::vector<TriangleContainer::SizeType> tList;
    tList = e1.getTriangles();
    BOOST_CHECK_MESSAGE(tList.size() == 0, "The assosiated triangles should be empty.");
}

BOOST_AUTO_TEST_CASE(associateAndGetTriangles)
{
    Node n1 {0.0, 0.0, 0.0};
    Node n2 {1.0, 1.0, 1.0};
    Edge e1 {n1, n2};
    e1.associateTriangle(3);
    e1.associateTriangle(6);
    e1.associateTriangle(2);
    e1.associateTriangle(4);

    std::vector<TriangleContainer::SizeType> tList;
    tList = e1.getTriangles();
    BOOST_CHECK_MESSAGE(tList.size() == 4, "There should be four triangle indeces associated with the edge.");
    BOOST_CHECK_MESSAGE(tList.at(0) == 2, "The first element should be 2 (sorted).");
    BOOST_CHECK_MESSAGE(tList.at(1) == 3, "The first element should be 3 (sorted).");
    BOOST_CHECK_MESSAGE(tList.at(2) == 4, "The first element should be 4 (sorted).");
    BOOST_CHECK_MESSAGE(tList.at(3) == 6, "The first element should be 6 (sorted).");

    e1.associateTriangle(6);
    e1.associateTriangle(3);
    e1.associateTriangle(1);
    tList = e1.getTriangles();
    BOOST_CHECK_MESSAGE(tList.size() == 5, "There should be five triangle indeces associated with the edge.");
    BOOST_CHECK_MESSAGE(tList.at(0) == 1, "The first element should be 1 (sorted + unique).");
    BOOST_CHECK_MESSAGE(tList.at(1) == 2, "The first element should be 2 (sorted + unique).");
    BOOST_CHECK_MESSAGE(tList.at(2) == 3, "The first element should be 3 (sorted + unique).");
    BOOST_CHECK_MESSAGE(tList.at(3) == 4, "The first element should be 4 (sorted + unique).");
    BOOST_CHECK_MESSAGE(tList.at(4) == 6, "The first element should be 6 (sorted + unique).");

}

BOOST_AUTO_TEST_SUITE_END()

//EOF

