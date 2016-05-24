#include "Node.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Node_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    Node p1;
    BOOST_CHECK_MESSAGE(p1.x() == 0.0, "The default value for point should be (0.0, 0.0, 0.0).");
    BOOST_CHECK_MESSAGE(p1.y() == 0.0, "The default value for point should be (0.0, 0.0, 0.0).");
    BOOST_CHECK_MESSAGE(p1.z() == 0.0, "The default value for point should be (0.0, 0.0, 0.0).");

    Node p2 {1.0, 2.0, 3.0};
    BOOST_CHECK_MESSAGE(p2.x() == 1.0, "The point was constructed with x = 1.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 2.0, "The point was constructed with y = 2.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 3.0, "The point was constructed with z = 3.0.");
}

BOOST_AUTO_TEST_CASE(setAndGetPoint)
{
    Node p1;
    p1.set(5.0, 15.0, 25.0);
    BOOST_CHECK_MESSAGE(p1.x() == 5.0 , "The point should have x = 5.0.");
    BOOST_CHECK_MESSAGE(p1.y() == 15.0, "The point should have y = 15.0.");
    BOOST_CHECK_MESSAGE(p1.z() == 25.0, "The point should have z = 25.0.");

}

BOOST_AUTO_TEST_CASE(pointEquality)
{
    Node p1;
    Node p2;

    BOOST_CHECK_MESSAGE(p1 == p2, "Points are expected to be equal.");

    p1.set(1.0, 2.0, 3.0);
    BOOST_CHECK_MESSAGE(not (p1 == p2), "Points are not expected to be equal.");
    BOOST_CHECK_MESSAGE(p1 != p2, "Points are not expected to be equal.");

    p2.set(1.0, 2.0, 3.0);
    BOOST_CHECK_MESSAGE(p1 == p2, "Points are expected to be equal.");
    BOOST_CHECK_MESSAGE(not (p1 != p2), "Points are expected to be equal.");
}


BOOST_AUTO_TEST_SUITE_END()

//EOF

