#include "NearFieldContainer.hpp"
#include "Node.hpp"
#include "NearFieldValue.hpp"
#include <iostream>

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(NearFieldContainer_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    NearFieldContainer fContainer;
    BOOST_CHECK_MESSAGE(fContainer.getFieldType() == NearFieldContainer::UNDEFINED, "Default field type is undefined.");
    BOOST_CHECK_MESSAGE(fContainer.size() == 0, "Container should be empty by default.");
}

BOOST_AUTO_TEST_CASE(setAndGetFieldType)
{
    NearFieldContainer fContainer;
    BOOST_CHECK_MESSAGE(fContainer.getFieldType() == NearFieldContainer::UNDEFINED, "Default field type is undefined.");

    fContainer.setFieldType(NearFieldContainer::ELECTICFIELD);
    BOOST_CHECK_MESSAGE(fContainer.getFieldType() == NearFieldContainer::ELECTICFIELD, "Should be electric field.");

    fContainer.setFieldType(NearFieldContainer::MAGNETICFIELD);
    BOOST_CHECK_MESSAGE(fContainer.getFieldType() == NearFieldContainer::MAGNETICFIELD, "Should be magnetic field.");
}

BOOST_AUTO_TEST_CASE(setPoint)
{
    NearFieldContainer fContainer;
    Node p1 {1.0, 2.0, 3.0};
    fContainer.setPoint(p1);
    BOOST_CHECK_MESSAGE(fContainer.size() == 1, "Should have a single point.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0).x() == p1.x(), "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0).y() == p1.y(), "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0).z() == p1.z(), "Incorrect value.");

    Node p2 {-1.0, -2.0, -3.0};
    fContainer.setPoint(p2);
    BOOST_CHECK_MESSAGE(fContainer.size() == 1, "Should have a single point.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0).x() == p2.x(), "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0).y() == p2.y(), "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0).z() == p2.z(), "Incorrect value.");
}

BOOST_AUTO_TEST_CASE(setPoints)
{
    NearFieldContainer fContainer;
    Node p1 {1.0, 2.0, 3.0};
    Node p2 {11.0, 4.0, 3.0};
    fContainer.setPoints(p1, p2, 11, 3, 1);
    BOOST_CHECK_MESSAGE(fContainer.size() == 33, "Should have 11x3x1 points.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(0) == p1, "Should be the first point.");
    BOOST_CHECK_MESSAGE(fContainer.getPointAt(32) == p2, "Should be the last point.");
}

BOOST_AUTO_TEST_CASE(setAndGetValues)
{
    NearFieldContainer fContainer;
    Node p1 {1.0, 2.0, 3.0};
    Node p2 {2.0, 3.0, 3.0};
    fContainer.setPoints(p1, p2, 2, 2, 1);
    NearFieldValue val1;
    NearFieldValue val2;
    NearFieldValue val3;
    NearFieldValue val4;
    val1.setX({1.0, -1.0});
    val2.setX({2.0, -2.0});
    val3.setX({3.0, -3.0});
    val4.setX({4.0, -4.0});
    fContainer.setValueAt(0, val1);
    fContainer.setValueAt(1, val2);
    fContainer.setValueAt(2, val3);
    fContainer.setValueAt(3, val4);
    BOOST_CHECK_MESSAGE(fContainer.getValueAt(0) == val1, "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getValueAt(1) == val2, "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getValueAt(2) == val3, "Incorrect value.");
    BOOST_CHECK_MESSAGE(fContainer.getValueAt(3) == val4, "Incorrect value.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

