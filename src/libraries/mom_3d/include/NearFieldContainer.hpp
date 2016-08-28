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
        void		setValueAt(unsigned index, NearFieldValue value);
        NearFieldValue 	getValueAt(unsigned index) const;
        unsigned	size() const;
        void		writeToEFEHFE(std::string fname);

    protected:
    private:
        unsigned	m_numPoints;
        FieldType	m_fieldType;
        NodeContainer 	m_nodes;
        std::vector<std::complex<double>> m_fieldValueXcomponent;
        std::vector<std::complex<double>> m_fieldValueYcomponent;
        std::vector<std::complex<double>> m_fieldValueZcomponent;

        double 		calcDelta(double start, double end, unsigned count);
};
