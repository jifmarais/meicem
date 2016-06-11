#include "Vector.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Vector_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    Vector p1;
    BOOST_CHECK_MESSAGE(p1.x() == 0.0, "The default value for point should be (0.0, 0.0, 0.0).");
    BOOST_CHECK_MESSAGE(p1.y() == 0.0, "The default value for point should be (0.0, 0.0, 0.0).");
    BOOST_CHECK_MESSAGE(p1.z() == 0.0, "The default value for point should be (0.0, 0.0, 0.0).");

    Vector p2 {1.0, 2.0, 3.0};
    BOOST_CHECK_MESSAGE(p2.x() == 1.0, "The point was constructed with x = 1.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 2.0, "The point was constructed with y = 2.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 3.0, "The point was constructed with z = 3.0.");
}

BOOST_AUTO_TEST_CASE(setAndGetPoint)
{
    Vector p1;
    p1.set(5.0, 15.0, 25.0);
    BOOST_CHECK_MESSAGE(p1.x() == 5.0 , "The point should have x = 5.0.");
    BOOST_CHECK_MESSAGE(p1.y() == 15.0, "The point should have y = 15.0.");
    BOOST_CHECK_MESSAGE(p1.z() == 25.0, "The point should have z = 25.0.");

}

BOOST_AUTO_TEST_CASE(pointEquality)
{
    Vector p1;
    Vector p2;

    BOOST_CHECK_MESSAGE(p1 == p2, "Points are expected to be equal.");

    p1.set(1.0, 2.0, 3.0);
    BOOST_CHECK_MESSAGE(not (p1 == p2), "Points are not expected to be equal.");
    BOOST_CHECK_MESSAGE(p1 != p2, "Points are not expected to be equal.");

    p2.set(1.0, 2.0, 3.0);
    BOOST_CHECK_MESSAGE(p1 == p2, "Points are expected to be equal.");
    BOOST_CHECK_MESSAGE(not (p1 != p2), "Points are expected to be equal.");
}

BOOST_AUTO_TEST_CASE(setEqaulVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);
    p2.set(3.0, 2.0, 1.0);

    p2 = p1;
    BOOST_CHECK_MESSAGE(p2.x() == 1.0 , "The point should have x = 1.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 2.0, "The point should have y = 2.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 3.0, "The point should have z = 3.0.");

    p2 = p2;
    BOOST_CHECK_MESSAGE(p2.x() == 1.0 , "The point should have x = 1.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 2.0, "The point should have y = 2.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 3.0, "The point should have z = 3.0.");
}

BOOST_AUTO_TEST_CASE(addEqualsVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);
    p2.set(3.0, 2.0, 1.0);

    p2 += p1;
    BOOST_CHECK_MESSAGE(p2.x() == 4.0 , "The point should have x = 4.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 4.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 4.0, "The point should have z = 4.0.");

    p2 += p1;
    BOOST_CHECK_MESSAGE(p2.x() == 5.0 , "The point should have x = 5.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 6.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 7.0, "The point should have z = 7.0.");
}

BOOST_AUTO_TEST_CASE(addVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);
    p2.set(3.0, 2.0, 1.0);

    p2 = p1 + p2;
    BOOST_CHECK_MESSAGE(p2.x() == 4.0 , "The point should have x = 4.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 4.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 4.0, "The point should have z = 4.0.");

    p2 = p1 + p2;
    BOOST_CHECK_MESSAGE(p2.x() == 5.0 , "The point should have x = 5.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 6.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 7.0, "The point should have z = 7.0.");
}

BOOST_AUTO_TEST_CASE(subtractEqualsVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);
    p2.set(3.0, 2.0, 1.0);

    p2 -= p1;
    BOOST_CHECK_MESSAGE(p2.x() == 2.0 , "The point should have x = 2.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 0.0, "The point should have y = 0.0.");
    BOOST_CHECK_MESSAGE(p2.z() == -2.0, "The point should have z = -2.0.");

    p2 -= p1;
    BOOST_CHECK_MESSAGE(p2.x() == 1.0 , "The point should have x = 5.0.");
    BOOST_CHECK_MESSAGE(p2.y() == -2.0, "The point should have y = -2.0.");
    BOOST_CHECK_MESSAGE(p2.z() == -5.0, "The point should have z = -5.0.");
}

BOOST_AUTO_TEST_CASE(subtractVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);
    p2.set(3.0, 2.0, 1.0);

    p2 = p1 - p2;
    BOOST_CHECK_MESSAGE(p2.x() == -2.0 , "The point should have x = -2.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 0.0, "The point should have y = 0.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 2.0, "The point should have z = 2.0.");

    p2 = p1 - p2;
    BOOST_CHECK_MESSAGE(p2.x() == 3.0 , "The point should have x = 5.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 2.0, "The point should have y = -2.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 1.0, "The point should have z = -5.0.");
}


BOOST_AUTO_TEST_CASE(multiplyEqualsVector)
{
    Vector p1;
    p1.set(1.0, 2.0, 3.0);

    p1 *= 2.0;
    BOOST_CHECK_MESSAGE(p1.x() == 2.0 , "The point should have x = 2.0.");
    BOOST_CHECK_MESSAGE(p1.y() == 4.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p1.z() == 6.0, "The point should have z = 6.0.");

    p1 *= -3.0;
    BOOST_CHECK_MESSAGE(p1.x() == -6.0 , "The point should have x = -6.0.");
    BOOST_CHECK_MESSAGE(p1.y() == -12.0, "The point should have y = -12.0.");
    BOOST_CHECK_MESSAGE(p1.z() == -18.0, "The point should have z = -18.0.");
}

BOOST_AUTO_TEST_CASE(multiplyVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);

    p2 = p1 * 2.0;
    BOOST_CHECK_MESSAGE(p2.x() == 2.0 , "The point should have x = 2.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 4.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 6.0, "The point should have z = 6.0.");

    p2 = p2 * -3.0;
    BOOST_CHECK_MESSAGE(p2.x() == -6.0 , "The point should have x = -6.0.");
    BOOST_CHECK_MESSAGE(p2.y() == -12.0, "The point should have y = -12.0.");
    BOOST_CHECK_MESSAGE(p2.z() == -18.0, "The point should have z = -18.0.");

    p2 = 2.0 * p2;
    BOOST_CHECK_MESSAGE(p2.x() == -12.0, "The point should have x = -12.0.");
    BOOST_CHECK_MESSAGE(p2.y() == -24.0, "The point should have y = -24.0.");
    BOOST_CHECK_MESSAGE(p2.z() == -36.0, "The point should have z = -36.0.");
}

BOOST_AUTO_TEST_CASE(divideEqualsVector)
{
    Vector p1;
    p1.set(1.0, 2.0, 3.0);

    p1 /= 1.0/2.0;
    BOOST_CHECK_MESSAGE(p1.x() == 2.0 , "The point should have x = 2.0.");
    BOOST_CHECK_MESSAGE(p1.y() == 4.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p1.z() == 6.0, "The point should have z = 6.0.");

    p1 /= -1.0/3.0;
    BOOST_CHECK_MESSAGE(p1.x() == -6.0 , "The point should have x = -6.0.");
    BOOST_CHECK_MESSAGE(p1.y() == -12.0, "The point should have y = -12.0.");
    BOOST_CHECK_MESSAGE(p1.z() == -18.0, "The point should have z = -18.0.");
}

BOOST_AUTO_TEST_CASE(divideVector)
{
    Vector p1;
    Vector p2;
    p1.set(1.0, 2.0, 3.0);

    p2 = p1 / (1.0/2.0);
    BOOST_CHECK_MESSAGE(p2.x() == 2.0 , "The point should have x = 2.0.");
    BOOST_CHECK_MESSAGE(p2.y() == 4.0, "The point should have y = 4.0.");
    BOOST_CHECK_MESSAGE(p2.z() == 6.0, "The point should have z = 6.0.");

    p2 = p2 / (-1.0/3.0);
    BOOST_CHECK_MESSAGE(p2.x() == -6.0 , "The point should have x = -6.0.");
    BOOST_CHECK_MESSAGE(p2.y() == -12.0, "The point should have y = -12.0.");
    BOOST_CHECK_MESSAGE(p2.z() == -18.0, "The point should have z = -18.0.");
}

BOOST_AUTO_TEST_CASE(magnitudeVector)
{
    Vector p1;
    p1.set(4.0, 0.0, 3.0);
    BOOST_CHECK_MESSAGE(p1.magnitude() == 5.0 , "The magnitude should be 5.0.");
    p1.set(4.0, 3.0, 0.0);
    BOOST_CHECK_MESSAGE(p1.magnitude() == 5.0 , "The magnitude should be 5.0.");
    p1.set(0.0, 0.0, 0.0);
    BOOST_CHECK_MESSAGE(p1.magnitude() == 0.0 , "The magnitude should be 0.0.");
    p1.set(0.0, -4.0, -3.0);
    BOOST_CHECK_MESSAGE(p1.magnitude() == 5.0 , "The magnitude should be 5.0.");
}

BOOST_AUTO_TEST_CASE(crossProductVector)
{
    Vector p1;
    Vector p2;
    Vector p3;

    p1.set(1.0, 0.0, 0.0);
    p2.set(0.0, 1.0, 0.0);
    p3 = Vector::cross(p1, p2);
    BOOST_CHECK_MESSAGE(p3.x() == 0.0, "The point should have x = 0.0.");
    BOOST_CHECK_MESSAGE(p3.y() == 0.0, "The point should have y = 0.0.");
    BOOST_CHECK_MESSAGE(p3.z() == 1.0, "The point should have z = 0.0.");

    p1.set(1.0, 0.0, 0.0);
    p2.set(0.0, 0.0, 1.0);
    p3 = Vector::cross(p1, p2);
    BOOST_CHECK_MESSAGE(p3.x() == 0.0, "The point should have x = 0.0.");
    BOOST_CHECK_MESSAGE(p3.y() == -1.0, "The point should have y = -1.0.");
    BOOST_CHECK_MESSAGE(p3.z() == 0.0, "The point should have z = 0.0.");
}

BOOST_AUTO_TEST_CASE(normVector)
{
    Vector p1;
    p1.set(4.0, 0.0, 3.0);
    BOOST_CHECK_MESSAGE(p1.norm().x() == 4.0/5.0 , "The magnitude should be 4.0/5.0.");
    BOOST_CHECK_MESSAGE(p1.norm().y() == 0.0 , "The magnitude should be 0.0.");
    BOOST_CHECK_MESSAGE(p1.norm().z() == 3.0/5.0 , "The magnitude should be 3.0/5.0.");
    p1.set(0.0, 0.0, 0.0);
    BOOST_CHECK_MESSAGE(p1.norm().x() == 0.0 , "The magnitude should be 0.0.");
    BOOST_CHECK_MESSAGE(p1.norm().y() == 0.0 , "The magnitude should be 0.0.");
    BOOST_CHECK_MESSAGE(p1.norm().z() == 0.0 , "The magnitude should be 0.0.");

    p1.set(-4.0, 3.0, 0.0);
    BOOST_CHECK_MESSAGE(p1.norm().x() == -4.0/5.0 , "The magnitude should be -4.0/5.0.");
    BOOST_CHECK_MESSAGE(p1.norm().y() == 3.0/5.0 , "The magnitude should be 30./5.0.");
    BOOST_CHECK_MESSAGE(p1.norm().z() == 0.0 , "The magnitude should be 0.0.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

