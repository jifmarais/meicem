#include "LabelContainer.hpp"
#include <algorithm>

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(LabelContainer_BasicTests)

BOOST_AUTO_TEST_CASE(addAndFindOneLabel)
{
    LabelContainer lContainer;
    lContainer.add("label1", 1);

    std::vector<LabelContainer::SizeType> iList;
    iList = lContainer.find("label1");
    BOOST_CHECK_MESSAGE(iList.size() == 1, "The list should only contain one element.");
    BOOST_CHECK_MESSAGE(iList[0] == 1, "The index should be 1.");
}

BOOST_AUTO_TEST_CASE(findMissingLabels)
{
    LabelContainer lContainer;

    std::vector<LabelContainer::SizeType> iList;
    iList = lContainer.find("label1");
    BOOST_CHECK_MESSAGE(iList.size() == 0, "The list should by empty.");

    lContainer.add("label1", 1);
    iList = lContainer.find("label2");
    BOOST_CHECK_MESSAGE(iList.size() == 0, "The list should by empty.");
}

BOOST_AUTO_TEST_CASE(addIndexMultipleTimes)
{
    LabelContainer lContainer;
    lContainer.add("label1", 1);
    lContainer.add("label1", 1);

    std::vector<LabelContainer::SizeType> iList;
    iList = lContainer.find("label1");
    BOOST_CHECK_MESSAGE(iList.size() == 1, "The list should only contain one element.");
    BOOST_CHECK_MESSAGE(iList[0] == 1, "The index should be 1.");
}

BOOST_AUTO_TEST_CASE(addIndexInMultipleLabel)
{
    LabelContainer lContainer;
    lContainer.add("label1", 1);
    lContainer.add("label2", 1);

    std::vector<LabelContainer::SizeType> iList;
    iList = lContainer.find("label1");
    BOOST_CHECK_MESSAGE(iList.size() == 1, "The list should only contain one element.");
    BOOST_CHECK_MESSAGE(iList[0] == 1, "The index should be 1.");

    iList = lContainer.find("label2");
    BOOST_CHECK_MESSAGE(iList.size() == 1, "The list should only contain one element.");
    BOOST_CHECK_MESSAGE(iList[0] == 1, "The index should be 1.");
}

BOOST_AUTO_TEST_CASE(findIndicesThatExist)
{
    LabelContainer lContainer;
    std::vector<LabelContainer::SizeType> iList;

    lContainer.add("label1", 1);
    lContainer.add("label1", 10);
    lContainer.add("label1", 100);

    lContainer.add("label2", 2);
    lContainer.add("label2", 20);
    lContainer.add("label2", 200);
    lContainer.add("label2", 2000);

    iList = lContainer.find("label1");
    BOOST_CHECK_MESSAGE(iList.size() == 3, "There should be 3 intries in the list.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 1) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 10) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 100) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 2) == iList.end(), "Index entry should not be found.");

    iList = lContainer.find("label2");
    BOOST_CHECK_MESSAGE(iList.size() == 4, "There should be 4 entries in the list.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 2) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 20) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 200) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 2000) != iList.end(), "Index entry should be found.");
    BOOST_CHECK_MESSAGE(std::find(iList.begin(), iList.end(), 1) == iList.end(), "Index entry should not be found.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

