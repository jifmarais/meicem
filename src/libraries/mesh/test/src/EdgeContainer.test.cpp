#include "Node.hpp"
#include "NodeContainer.hpp"
#include "Triangle.hpp"
#include "EdgeContainer.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(EdgeContainer_BasicTests)

BOOST_AUTO_TEST_CASE(constructorEmpty)
{
    Node n1;
    Node n2;
    Node n3;
    Triangle t;
    NodeContainer nContainer;
    TriangleContainer tContainer(nContainer);

    EdgeContainer eContainer(tContainer);
    eContainer.buildNonboundaryEdgeList();
    BOOST_CHECK_MESSAGE(eContainer.size() == 0, "No triangles = no edges.");
}

BOOST_AUTO_TEST_CASE(constructorOnlyBoundaryEdges)
{
    Node n1;
    Node n2;
    Node n3;
    Triangle t;
    NodeContainer nContainer;
    TriangleContainer tContainer(nContainer);

    n1.set(0.0, 0.0, 0.0);
    n2.set(1.0, 0.0, 0.0);
    n3.set(0.0, 1.0, 0.0);
    t.set(n1, n2, n3);
    tContainer.add(t);

    EdgeContainer eContainer(tContainer);
    eContainer.buildNonboundaryEdgeList();
    BOOST_CHECK_MESSAGE(eContainer.size() == 0, "One triangle = no edges.");

    n1.set(0.0, 0.0, 10.0);
    n2.set(10.0, 0.0, 0.0);
    n3.set(0.0, 10.0, 0.0);
    t.set(n1, n2, n3);
    tContainer.add(t);

    eContainer.buildNonboundaryEdgeList();
    BOOST_CHECK_MESSAGE(eContainer.size() == 0, "Two non-touching tiangles = no edges.");
}

BOOST_AUTO_TEST_CASE(constructor)
{
    Node n1;
    Node n2;
    Node n3;
    Triangle t;
    NodeContainer nContainer;
    TriangleContainer tContainer(nContainer);

    n1.set(0.0, 0.0, 0.0);
    n2.set(1.0, 0.0, 0.0);
    n3.set(0.0, 1.0, 0.0);
    t.set(n1, n2, n3);
    TriangleContainer::SizeType i1;
    i1 = tContainer.add(t);

    n1.set(0.0, 0.0, 0.0);
    n2.set(-1.0, 0.0, 0.0);
    n3.set(0.0, 1.0, 0.0);
    t.set(n1, n2, n3);
    TriangleContainer::SizeType i2;
    i2 = tContainer.add(t);

    n1.set(0.0, 0.0, 0.0);
    n2.set(1.0, 0.0, 0.0);
    n3.set(0.0, -1.0, 0.0);
    t.set(n1, n2, n3);
    TriangleContainer::SizeType i3;
    i3 = tContainer.add(t);

    EdgeContainer eContainer(tContainer);
    eContainer.buildNonboundaryEdgeList();
    BOOST_CHECK_MESSAGE(eContainer.size() == 2, "Should have two non-boundary edges.");
    BOOST_CHECK_MESSAGE(eContainer.at(0).getTriangles().size() == 2, "Should have two triangles associated.");
    BOOST_CHECK_MESSAGE(eContainer.at(0).getTriangles().at(0) == i1, "Should have the first triangle associated.");
    BOOST_CHECK_MESSAGE(eContainer.at(0).getTriangles().at(1) == i2, "Should have the second triangle associated.");
    BOOST_CHECK_MESSAGE(eContainer.at(1).getTriangles().size() == 2, "Should have two triangles associated.");
    BOOST_CHECK_MESSAGE(eContainer.at(1).getTriangles().at(0) == i1, "Should have the first triangle associated.");
    BOOST_CHECK_MESSAGE(eContainer.at(1).getTriangles().at(1) == i3, "Should have the second triangle associated.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

