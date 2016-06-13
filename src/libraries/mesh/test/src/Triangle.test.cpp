#include "Triangle.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Triangle_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    Node p0 {0.0, 0.0, 0.0};
    Triangle t1;
    BOOST_CHECK_MESSAGE(t1.n1() == p0, "The default indices for nodes should be the origin.");
    BOOST_CHECK_MESSAGE(t1.n2() == p0, "The default indices for nodes should be the origin.");
    BOOST_CHECK_MESSAGE(t1.n3() == p0, "The default indices for nodes should be the origin.");

    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    Triangle t2 {p1, p2, p3};
    BOOST_CHECK_MESSAGE(t2.n1() == p1, "The triangle was constructed with n1 = p1.");
    BOOST_CHECK_MESSAGE(t2.n2() == p2, "The triangle was constructed with n2 = p2.");
    BOOST_CHECK_MESSAGE(t2.n3() == p3, "The triangle was constructed with n3 = p3.");
}

BOOST_AUTO_TEST_CASE(setAndGetTriangle)
{
    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    Triangle t2;
    t2.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t2.n1() == p1, "The triangle should have n1 = p1.");
    BOOST_CHECK_MESSAGE(t2.n2() == p2, "The triangle should have n2 = p2.");
    BOOST_CHECK_MESSAGE(t2.n3() == p3, "The triangle should have n3 = p3.");
}

BOOST_AUTO_TEST_CASE(triangleEquality)
{
    Triangle t1;
    Triangle t2;
    BOOST_CHECK_MESSAGE(t1 == t2, "Triangles are expected to be equal.");

    Node p1 {1.0, 2.0, 3.0};
    Node p2 {4.0, 5.0, 6.0};
    Node p3 {7.0, 8.0, 9.0};
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(not (t1 == t2), "Triangles are not expected to be equal.");
    BOOST_CHECK_MESSAGE(t1 != t2, "Triangles are not expected to be equal.");

    t2.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1 == t2, "Triangles are expected to be equal.");
    BOOST_CHECK_MESSAGE(not (t1 != t2), "Triangles are expected to be equal.");
}

BOOST_AUTO_TEST_CASE(triangleArea)
{
    Triangle t1;
    Node p1 {0.0, 0.0, 0.0};
    Node p2 {0.0, 0.0, 0.0};
    Node p3 {0.0, 0.0, 0.0};
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1.area() == 0.0, "Traingle area should be 0.0.");

    p1.set(0.0, 0.0, 0.0);
    p2.set(1.0, 0.0, 0.0);
    p3.set(1.0, 1.0, 0.0);
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1.area() == 0.5, "Traingle area should be 0.5.");

    p1.set(0.0, 0.0, 0.0);
    p2.set(2.0, 0.0, 0.0);
    p3.set(2.0, 3.0, 4.0);
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1.area() == 5.0, "Traingle area should be 5.");

    p3.set(2.0, 0.0, 0.0);
    p2.set(-2.0, 0.0, 0.0);
    p1.set(-2.0, -3.0, -4.0);
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1.area() == 10.0, "Traingle area should be 10.");
}

BOOST_AUTO_TEST_CASE(triangleCentre)
{
    Triangle t1;
    Node p1 {0.0, 0.0, 0.0};
    Node p2 {0.0, 0.0, 0.0};
    Node p3 {0.0, 0.0, 0.0};
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1.centre().x() == 0.0, "Traingle centre should be 0.0.");
    BOOST_CHECK_MESSAGE(t1.centre().y() == 0.0, "Traingle centre should be 0.0.");
    BOOST_CHECK_MESSAGE(t1.centre().z() == 0.0, "Traingle centre should be 0.0.");

    p1.set(0.0, 0.0, 1.0);
    p2.set(1.0, 0.0, 1.0);
    p3.set(1.0, 1.0, 1.0);
    t1.set(p1, p2, p3);
    BOOST_CHECK_MESSAGE(t1.centre().x() == 2.0/3.0, "Traingle centre should be 2.0/3.0.");
    BOOST_CHECK_MESSAGE(t1.centre().y() == 1.0/3.0, "Traingle centre should be 1.0/3.0.");
    BOOST_CHECK_MESSAGE(t1.centre().z() == 1.0, "Traingle centre should be 1.0.");
}

BOOST_AUTO_TEST_CASE(triangleNormal)
{
    Triangle t1;
    Node p1 {0.0, 0.0, 0.0};
    Node p2 {2.0, 0.0, 0.0};
    Node p3 {0.0, 2.0, 0.0};
    Node answer;

    t1.set(p1, p2, p3);
    answer.set(0.0, 0.0, 1.0);
    BOOST_CHECK_MESSAGE(t1.normal() == answer, "Normal vector not correct.");

    p1.set(0.0, 0.0, 2.0);
    t1.set(p1, p2, p3);
    answer.set(1.0, 1.0, 1.0);
    BOOST_CHECK_MESSAGE(t1.normal() == answer.norm(), "Normal vector not correct.");
}

BOOST_AUTO_TEST_CASE(toSimplex)
{
    // TODO: I'm not happy wiht this implemenation. Too many open issues.
    Triangle t1;
    Node p1 {0.0, 0.0, 2.0};
    Node p2 {2.0, 0.0, 0.0};
    Node p3 {0.0, 2.0, 0.0};
    t1.set(p1, p2, p3);

    Node testPoint {0.0, 0.0, 2.0};
    Node answer {1.0, 0.0, 0.0};
    BOOST_CHECK_MESSAGE(t1.toSimplex(testPoint) == answer, "Conversion to simplex failed.");

    testPoint.set(2.0, 0.0, 0.0);
    answer.set(0.0, 1.0, 0.0);
    BOOST_CHECK_MESSAGE(t1.toSimplex(testPoint) == answer, "Conversion to simplex failed.");

    testPoint.set(0.0, 2.0, 0.0);
    answer.set(0.0, 0.0, 1.0);
    BOOST_CHECK_MESSAGE(t1.toSimplex(testPoint) == answer, "Conversion to simplex failed.");

    testPoint.set(2.0/3.0, 2.0/3.0, 2.0/3.0);
    answer.set(1.0/3.0, 1.0/3.0, 1.0/3.0);
    BOOST_CHECK_MESSAGE(t1.toSimplex(testPoint) == answer, "Conversion to simplex failed.");

    // What about when the point is not in the plane of the triangle?
}

BOOST_AUTO_TEST_CASE(fromSimplex)
{
    Triangle t1;
    Node p1 {0.0, 0.0, 2.0};
    Node p2 {2.0, 0.0, 0.0};
    Node p3 {0.0, 2.0, 0.0};
    t1.set(p1, p2, p3);

    Node testPoint {0.0, 0.0, 1.0};
    Node answer {0.0, 2.0, 0.0};
    BOOST_CHECK_MESSAGE(t1.fromSimplex(testPoint) == answer, "Conversion from simplex failed.");

    testPoint.set(1.0, 0.0, 0.0);
    answer.set(0.0, 0.0, 2.0);
    BOOST_CHECK_MESSAGE(t1.fromSimplex(testPoint) == answer, "Conversion from simplex failed.");

    testPoint.set(0.0, 1.0, 0.0);
    answer.set(2.0, 0.0, 0.0);
    BOOST_CHECK_MESSAGE(t1.fromSimplex(testPoint) == answer, "Conversion from simplex failed.");

    testPoint.set(1.0/3.0, 1.0/3.0, 1.0/3.0);
    BOOST_CHECK_MESSAGE(t1.fromSimplex(testPoint) == t1.centre(), "Should be the triangle centre.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

