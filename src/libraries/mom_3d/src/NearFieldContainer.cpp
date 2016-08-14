#include "NearFieldContainer.hpp"
#include "assert.h"
#include "fstream"

NearFieldContainer::NearFieldContainer()
{
    //ctor
    m_numPoints = 0;
    m_fieldType = FieldType::UNDEFINED;
}

NearFieldContainer::~NearFieldContainer()
{
    //dtor
}

void NearFieldContainer::setFieldType(NearFieldContainer::FieldType field)
{
    m_fieldType = field;
}

double NearFieldContainer::calcDelta(double start, double end, unsigned count)
{
    assert(count > 0);

    double dx;
    if (count == 1)
    {
        dx = 0.0;
    }
    else
    {
        dx = (end - start) / (double)(count-1);
    }
    return dx;
}

void NearFieldContainer::setPoints(Node start,
                                   Node end,
                                   unsigned xCount,
                                   unsigned yCount,
                                   unsigned zCount)
{
    assert(xCount > 0);
    assert(yCount > 0);
    assert(zCount > 0);

    m_fieldValueXcomponent.clear();
    m_fieldValueYcomponent.clear();
    m_fieldValueZcomponent.clear();
    m_nodes.clear();

    double dx = calcDelta(start.x(), end.x(), xCount);
    double dy = calcDelta(start.y(), end.y(), yCount);
    double dz = calcDelta(start.y(), end.y(), zCount);

    Node p = start;
    for (auto xIndex = 0 ; xIndex < (int)xCount; ++xIndex)
    {
        double xCoord = start.x() + (double)xIndex*(dx);
        for (auto yIndex = 0 ; yIndex < (int)yCount; ++yIndex)
        {
            double yCoord = start.y() + (double)yIndex*(dy);
            for (auto zIndex = 0 ; zIndex < (int)zCount; ++zIndex)
            {
                double zCoord = start.z() + (double)zIndex*(dz);

                p.set(xCoord, yCoord, zCoord);
                m_nodes.add(p);
                m_fieldValueXcomponent.push_back({0.0, 0.0});
                m_fieldValueYcomponent.push_back({0.0, 0.0});
                m_fieldValueZcomponent.push_back({0.0, 0.0});
            }
        }
    }
}

void NearFieldContainer::setPoint(Node point)
{
    m_fieldValueXcomponent.clear();
    m_fieldValueYcomponent.clear();
    m_fieldValueZcomponent.clear();
    m_nodes.clear();

    m_nodes.add(point);
    m_fieldValueXcomponent.push_back({0.0, 0.0});
    m_fieldValueYcomponent.push_back({0.0, 0.0});
    m_fieldValueZcomponent.push_back({0.0, 0.0});
}

unsigned NearFieldContainer::size() const
{
    assert(m_nodes.size() == m_fieldValueXcomponent.size());
    assert(m_nodes.size() == m_fieldValueYcomponent.size());
    assert(m_nodes.size() == m_fieldValueZcomponent.size());

    return m_nodes.size();
}

NearFieldContainer::FieldType NearFieldContainer::getFieldType() const
{
    return m_fieldType;
}

void NearFieldContainer::setValueAt(unsigned index, NearFieldValue value)
{
    assert(index < size());

    m_fieldValueXcomponent.at(index) = value.getX();
    m_fieldValueYcomponent.at(index) = value.getY();
    m_fieldValueZcomponent.at(index) = value.getZ();
}

Node NearFieldContainer::getPointAt(unsigned index) const
{
    assert(index < size());
    return m_nodes.at(index);
}

NearFieldValue NearFieldContainer::getValueAt(unsigned index) const
{
    assert(index < size());

    NearFieldValue result;
    result.setX(m_fieldValueXcomponent.at(index));
    result.setY(m_fieldValueYcomponent.at(index));
    result.setZ(m_fieldValueZcomponent.at(index));
    return result;
}

void NearFieldContainer::writeToEFEHFE(std::string fname)
{
    assert(m_fieldType != NearFieldContainer::UNDEFINED);
    assert(fname != "");

    std::string ext;
    if ( NearFieldContainer::ELECTICFIELD == m_fieldType )
    {
        ext = "efe";
    }
    else
    {
        ext = "hfe";
    }
    std::ofstream file (fname + '.' + ext);

    if (file.is_open())
    {
        file << "JIF" << std::endl;
        file.close();
    }

}
