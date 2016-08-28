#include "NearFieldContainer.hpp"
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include "assert.h"
#include "fstream"

NearFieldContainer::NearFieldContainer()
{
    //ctor
    m_numPointsX = 0;
    m_numPointsY = 0;
    m_numPointsZ = 0;
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

double NearFieldContainer::calcDelta(double start, double end, unsigned count) const
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

unsigned NearFieldContainer::localToGlobalIndex(unsigned xIndex, unsigned yIndex, unsigned zIndex) const
{
//    return xIndex*(m_numPointsY*m_numPointsZ) + yIndex*m_numPointsZ + zIndex;
    return xIndex*(m_numPointsY*m_numPointsZ) + yIndex*m_numPointsZ + zIndex;
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

    m_numPointsX = xCount;
    m_numPointsY = yCount;
    m_numPointsZ = zCount;

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

void NearFieldContainer::setValueAt(unsigned xIndex, unsigned yIndex, unsigned zIndex , NearFieldValue value)
{
    assert(xIndex < m_numPointsX);
    assert(yIndex < m_numPointsY);
    assert(zIndex < m_numPointsZ);

    setValueAt(localToGlobalIndex(xIndex, yIndex, zIndex), value);
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

Node NearFieldContainer::getPointAt(unsigned xIndex, unsigned yIndex, unsigned zIndex) const
{
    assert(xIndex < m_numPointsX);
    assert(yIndex < m_numPointsY);
    assert(zIndex < m_numPointsZ);

    return getPointAt(localToGlobalIndex(xIndex, yIndex, zIndex));
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

NearFieldValue NearFieldContainer::getValueAt(unsigned xIndex, unsigned yIndex, unsigned zIndex) const
{
    assert(xIndex < m_numPointsX);
    assert(yIndex < m_numPointsY);
    assert(zIndex < m_numPointsZ);

    return getValueAt(localToGlobalIndex(xIndex, yIndex, zIndex));
}

void NearFieldContainer::writeToEFEHFE(std::string fname) const
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
        if (m_fieldType == NearFieldContainer::ELECTICFIELD)
        {
            file << "##File Type: Electric near field" << std::endl;
        }
        else
        {
            file << "##File Type: Magnetic near field" << std::endl;
        }
        file << "##File Format: 4" << std::endl;
        file << "** File exported by MEICEM" << std::endl;
        file << std::endl;
//        file << "#Configuration Name: StandardConfiguration1" << std::endl;
//        file << "#Request Name: NearField11" << std::endl;
        file << "#Frequency:   1.00000000E+009" << std::endl;
        file << "#Coordinate System: Cartesian" << std::endl;
        file << "#No. of X Samples: " << m_numPointsX << std::endl;
        file << "#No. of Y Samples: " << m_numPointsY << std::endl;
        file << "#No. of Z Samples: " << m_numPointsZ << std::endl;
        if (m_fieldType == NearFieldContainer::ELECTICFIELD)
        {
            file << "#Result Type: Electric Field Values" << std::endl;
        }
        else
        {
            file << "#Result Type: Magnetic Field Values" << std::endl;
        }
        file << "#No. of Header Lines: 1" << std::endl;
        file << "#         \"X\"       ";
        file << "        \"Y\"         ";
        file << "        \"Z\"         ";
        file << "      \"Re(Ex)\"      ";
        file << "      \"Im(Ex)\"      ";
        file << "      \"Re(Ey)\"      ";
        file << "      \"Im(Ey)\"      ";
        file << "      \"Re(Ez)\"      ";
        file << "      \"Im(Ez)\" ";
        file << std::endl;

        for (unsigned zIndex=0; zIndex < m_numPointsZ ; ++zIndex)
        {
            for (unsigned yIndex=0; yIndex < m_numPointsY ; ++yIndex)
            {
                for (unsigned xIndex=0; xIndex < m_numPointsX ; ++xIndex)
                {
                    file << std::fixed << std::setprecision(9) ;
                    file << getPointAt(xIndex, yIndex, zIndex).x() << "  " ;
                    file << getPointAt(xIndex, yIndex, zIndex).y() << "  " ;
                    file << getPointAt(xIndex, yIndex, zIndex).z() << "  " ;
                    file << getValueAt(xIndex, yIndex, zIndex).getX().real() << "  " ;
                    file << getValueAt(xIndex, yIndex, zIndex).getX().imag() << "  " ;
                    file << getValueAt(xIndex, yIndex, zIndex).getY().real() << "  " ;
                    file << getValueAt(xIndex, yIndex, zIndex).getY().imag() << "  " ;
                    file << getValueAt(xIndex, yIndex, zIndex).getZ().real() << "  " ;
                    file << getValueAt(xIndex, yIndex, zIndex).getZ().imag() ;
                    file << std::endl;
                }
            }
        }

        file << std::endl;
        file.close();
    }

}
