#pragma once
#include "NearFieldValue.hpp"
#include "NodeContainer.hpp"
#include "ComplexMatrix.hpp"

class NearFieldContainer
{
    public:
        enum	FieldType
                    {
                            UNDEFINED,
                            ELECTICFIELD,
                            MAGNETICFIELD
                    };

        NearFieldContainer();
        virtual         ~NearFieldContainer();
        void 		setFieldType(FieldType field);
        FieldType	getFieldType() const;
        void 		setPoints(Node   start,
                                  Node   end,
                                  unsigned xCount,
                                  unsigned yCount,
                                  unsigned zCount);
        void 		setPoint(Node point);
        Node		getPointAt(unsigned index) const;
        Node		getPointAt(unsigned xIndex, unsigned yIndex, unsigned zIndex) const;
        void		setValueAt(unsigned index, NearFieldValue value);
        void 		setValueAt(unsigned xIndex, unsigned yIndex, unsigned zIndex, NearFieldValue value);
        NearFieldValue 	getValueAt(unsigned index) const;
        NearFieldValue 	getValueAt(unsigned xIndex, unsigned yIndex, unsigned zIndex) const;
        unsigned	size() const;
        void		writeToEFEHFE(std::string fname) const;

protected:
    private:
        FieldType	m_fieldType;
        NodeContainer 	m_nodes;
        unsigned	m_numPointsX;
        unsigned	m_numPointsY;
        unsigned	m_numPointsZ;
        std::vector<std::complex<double>> m_fieldValueXcomponent;
        std::vector<std::complex<double>> m_fieldValueYcomponent;
        std::vector<std::complex<double>> m_fieldValueZcomponent;

        double 		calcDelta(double start, double end, unsigned count) const;
        unsigned	localToGlobalIndex(unsigned xIndex, unsigned yIndex, unsigned zIndex) const;
};
