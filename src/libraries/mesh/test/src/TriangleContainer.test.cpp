#include "Node.hpp"
#include "NodeContainer.hpp"
#include "Triangle.hpp"
#include "TriangleContainer.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(TriangleContainer_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    BOOST_CHECK_MESSAGE(tContainer.size() == 0, "The container should start empty.");
}

BOOST_AUTO_TEST_CASE(addTriangles)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);

    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    Triangle t1 {p1, p2, p3};
    NodeContainer::SizeType  i1 = tContainer.add(t1);
    BOOST_CHECK_MESSAGE(tContainer.size() == 1, "The container should have 1 triangle.");

    NodeContainer::SizeType  i2 = tContainer.add(t1);
    BOOST_CHECK_MESSAGE(tContainer.size() == 2, "The container should have 2 triangles.");

    //Currently duplicate triangles are allowed
    BOOST_CHECK_MESSAGE(i1 != i2, "Indices should differ even if the points are the same.");
}

BOOST_AUTO_TEST_CASE(triangleAtIndex)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);

    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    Triangle t1 {p1, p2, p3};
    NodeContainer::SizeType  i1 = tContainer.add(t1);

    p1.set(11.0, 12.0, 13.0);
    p2.set(14.0, 15.0, 16.0);
    p3.set(17.0, 18.0, 19.0);
    t1.set(p1, p2, p3);
    NodeContainer::SizeType  i2 = tContainer.add(t1);

    p1.set(21.0, 22.0, 23.0);
    p2.set(24.0, 25.0, 26.0);
    p3.set(27.0, 28.0, 29.0);
    t1.set(p1, p2, p3);
    NodeContainer::SizeType  i3 = tContainer.add(t1);

    t1 = tContainer.at(i2);
    BOOST_CHECK_MESSAGE(t1[0] == t1.n1()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1[1] == t1.n2()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1[2] == t1.n3()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1.n1().x() == 11.0, "The coordinate should be 11.");
    BOOST_CHECK_MESSAGE(t1.n2().y() == 15.0, "The coordinate should be 15.");
    BOOST_CHECK_MESSAGE(t1.n3().z() == 19.0, "The coordinate should be 19.");

    t1 = tContainer.at(i1);
    BOOST_CHECK_MESSAGE(t1[0] == t1.n1()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1[1] == t1.n2()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1[2] == t1.n3()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1.n1().y() == 2.0, "The coordinate should be 2.");
    BOOST_CHECK_MESSAGE(t1.n2().z() == 6.0, "The coordinate should be 6.");
    BOOST_CHECK_MESSAGE(t1.n3().x() == 7.0, "The coordinate should be 7.");

    t1 = tContainer.at(i3);
    BOOST_CHECK_MESSAGE(t1[0] == t1.n1()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1[1] == t1.n2()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1[2] == t1.n3()   , "Two ways to access the same node should be equal.");
    BOOST_CHECK_MESSAGE(t1.n1().z() == 23.0, "The coordinate should be 23.");
    BOOST_CHECK_MESSAGE(t1.n2().x() == 24.0, "The coordinate should be 24.");
    BOOST_CHECK_MESSAGE(t1.n3().y() == 28.0, "The coordinate should be 28.");

}

BOOST_AUTO_TEST_CASE(findPointByCoordinates)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);

    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    Triangle t1 {p1, p2, p3};
    NodeContainer::SizeType  i1 = tContainer.add(t1);

    p1.set(11.0, 12.0, 13.0);
    p2.set(14.0, 15.0, 16.0);
    p3.set(17.0, 18.0, 19.0);
    Triangle t2 {p1, p2, p3};
    NodeContainer::SizeType  i2 = tContainer.add(t2);

    p1.set(21.0, 22.0, 23.0);
    p2.set(24.0, 25.0, 26.0);
    p3.set(27.0, 28.0, 29.0);
    Triangle t3 {p1, p2, p3};
    NodeContainer::SizeType  i3 = tContainer.add(t3);

    BOOST_CHECK_MESSAGE(tContainer.find(t2) == i2, "Incorrect index found.");
    BOOST_CHECK_MESSAGE(tContainer.find(t1) == i1, "Incorrect index found.");
    BOOST_CHECK_MESSAGE(tContainer.find(t3) == i3, "Incorrect index found.");

    p2.set(124.0, 125.0, 26.0);
    Triangle t4 {p1, p2, p3};
    BOOST_CHECK_MESSAGE(tContainer.find(t4) == TriangleContainer::invalidIndex, "Triangles should not have been found.");

}

BOOST_AUTO_TEST_CASE(hasCommonNode)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);

    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    Triangle t1 {p1, p2, p3};
    NodeContainer::SizeType  i1 = tContainer.add(t1);

    p1.set(11.0, 12.0, 13.0);
    p2.set(14.0, 15.0, 16.0);
    p3.set(17.0, 18.0, 19.0);
    t1.set(p1, p2, p3);
    NodeContainer::SizeType  i2 = tContainer.add(t1);

    p1.set(1.0, 2.0, 3.0);
    t1.set(p1, p2, p3);
    NodeContainer::SizeType  i3 = tContainer.add(t1);

    BOOST_CHECK_MESSAGE(tContainer.hasCommonNode(i1, i2) == false, "Triangles don't share a node.");
    BOOST_CHECK_MESSAGE(tContainer.hasCommonNode(i1, i3) == true, "Triangles do share a node.");
    BOOST_CHECK_MESSAGE(tContainer.hasCommonNode(i2, i3) == true, "Triangles share multiple nodes.");

}

BOOST_AUTO_TEST_SUITE_END()

//EOF

