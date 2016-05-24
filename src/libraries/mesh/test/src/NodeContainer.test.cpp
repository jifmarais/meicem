#include "Node.hpp"
#include "NodeContainer.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(NodeContainer_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    NodeContainer L1;
    BOOST_CHECK_EQUAL(L1.size(), 0);
}

BOOST_AUTO_TEST_CASE(adds)
{
    NodeContainer L1;
    Node p1;
    NodeContainer::SizeType ii;
    NodeContainer::SizeType count {50};

    // Add a bunch of points
    for ( ii=0; ii < count ; ++ii )
    {
        p1.set(1.0*ii, 2.0*ii, 3.0*ii);
        L1.add(p1);
    }
    BOOST_CHECK_EQUAL(L1.size(), count);

    // Points should not be added multiple times 
    for ( ii=0; ii < count ; ++ii )
    {
        p1.set(1.0*ii, 2.0*ii, 3.0*ii);
        L1.add(p1);
    }
    BOOST_CHECK_EQUAL(L1.size(), count);
}

BOOST_AUTO_TEST_CASE(at)
{
    NodeContainer L1;
    Node p1;
    Node p2;
    NodeContainer::SizeType index;
    for ( NodeContainer::SizeType ii=0; ii <= 4 ; ++ii )
    {
        p1.set(1.0*ii, 2.0*ii, 3.0*ii);
        index = L1.add(p1);

        p2 = L1.at(index);
        BOOST_CHECK_MESSAGE(p2 == p1, "Points are expected to be equal.");
    }

    p1.set(1.0*2, 2.0*2, 3.0*2);
    p2 = L1.at(2);
    BOOST_CHECK_MESSAGE(p2 == p1, "Points are expected to be equal.");
}

BOOST_AUTO_TEST_CASE(findPoint)
{
    NodeContainer L1;
    Node p1;
    NodeContainer::SizeType index;
    for ( NodeContainer::SizeType ii=0; ii <= 4 ; ++ii )
    {
        p1.set(1.0*ii, 2.0*ii, 3.0*ii);
        index = L1.add(p1);

        BOOST_CHECK_MESSAGE(L1.find(p1) == index, "Incorrect index found.");
    }

    p1.set(1.0*2, 2.0*2, 3.0*2);
    BOOST_CHECK_MESSAGE(L1.find(p1) == 2, "Expecting to find index 2.");

    p1.set(1.0, 1.0, 1.0);
    BOOST_CHECK_MESSAGE(L1.find(p1) == NodeContainer::invalidIndex, "Not expecting a valid index to be returned." );
}

BOOST_AUTO_TEST_CASE(SetgetTolerance)
{
    NodeContainer L1;

    BOOST_CHECK_MESSAGE(L1.getTolerance() == 1e-6, "The default tolerance is not 1e-6.");

    L1.setTolerance(1e-5);
    BOOST_CHECK_MESSAGE(L1.getTolerance() == 1e-5, "The tolerance value is not what it was set to.");
}

BOOST_AUTO_TEST_CASE(ToleranceLimitsadd)
{
    NodeContainer L1;
    NodeContainer::SizeType index;
    NodeContainer::SizeType indexReference;

    double tolerance = 1e-2;
    L1.setTolerance(tolerance);

    Node p1 {1.0, 2.0, 3.0};
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 1);
    indexReference = index;

    p1.set(1.0 + 0.99*tolerance, 2.0, 3.0);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 1);
    BOOST_CHECK_EQUAL(indexReference, index);

    p1.set(1.0, 2.0 + 0.99*tolerance, 3.0);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 1);
    BOOST_CHECK_EQUAL(indexReference, index);

    p1.set(1.0, 2.0, 3.0 + 0.99*tolerance);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 1);
    BOOST_CHECK_EQUAL(indexReference, index);

    p1.set(1.0 + 1.01*tolerance, 2.0, 3.0);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 2);
    BOOST_CHECK_NE(indexReference, index);

    p1.set(1.0, 2.0 + 1.01*tolerance, 3.0);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 3);
    BOOST_CHECK_NE(indexReference, index);

    p1.set(1.0, 2.0, 3.0 + 1.01*tolerance);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 4);
    BOOST_CHECK_NE(indexReference, index);

}

BOOST_AUTO_TEST_CASE(ToleranceLimitsfindPoint)
{
    NodeContainer L1;
    NodeContainer::SizeType index;
    NodeContainer::SizeType indexReference;

    double tolerance = 1e-6;
    L1.setTolerance(tolerance);

    Node p1;
    p1.set(1.0, 2.0, 3.0);
    index = L1.add(p1);
    BOOST_CHECK_EQUAL(L1.size(), 1);
    indexReference = index;

    p1.set(1.0 + 0.99*tolerance, 2.0, 3.0);
    BOOST_CHECK_MESSAGE(L1.find(p1) == indexReference, "Incorrect index returned.");

    p1.set(1.0, 2.0 + 0.99*tolerance, 3.0);
    BOOST_CHECK_MESSAGE(L1.find(p1) == indexReference, "Incorrect index returned.");

    p1.set(1.0, 2.0, 3.0 + 0.99*tolerance);
    BOOST_CHECK_MESSAGE(L1.find(p1) == indexReference, "Incorrect index returned.");

    p1.set(1.0 + 1.01*tolerance, 2.0, 3.0);
    BOOST_CHECK_MESSAGE(L1.find(p1) == NodeContainer::invalidIndex, "Not expecting a valid index to be returned.");

    p1.set(1.0, 2.0 + 1.01*tolerance, 3.0);
    BOOST_CHECK_MESSAGE(L1.find(p1) == NodeContainer::invalidIndex, "Not expecting a valid index to be returned.");

    p1.set(1.0, 2.0, 3.0 + 1.01*tolerance);
    BOOST_CHECK_MESSAGE(L1.find(p1) == NodeContainer::invalidIndex, "Not expecting a valid index to be returned.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

