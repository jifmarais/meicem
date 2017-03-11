#include "NearFieldValue.hpp"
#include "assert.h"

NearFieldValue::NearFieldValue()
{
    //ctor
    m_field(0) = {0.0, 0.0};
    m_field(1) = {0.0, 0.0};
    m_field(2) = {0.0, 0.0};
}

NearFieldValue::~NearFieldValue()
{
    //dtor
}

void NearFieldValue::set(std::complex<double> x, std::complex<double> y, std::complex<double> z)
{
    m_field(0) = x;
    m_field(1) = y;
    m_field(2) = z;
}

void NearFieldValue::setX(std::complex<double> val)
{
    m_field(0) = val;
}

void NearFieldValue::setY(std::complex<double> val)
{
    m_field(1) = val;
}

void NearFieldValue::setZ(std::complex<double> val)
{
    m_field(2) = val;
}

NearFieldValue::cvec NearFieldValue::getVec() const
{
    return m_field;
}

void NearFieldValue::setVec(NearFieldValue::cvec value)
{
    m_field = value;
}


std::complex<double> NearFieldValue::getX() const
{
    return m_field(0);
}

std::complex<double> NearFieldValue::getY() const
{
    return m_field(1);
}

std::complex<double> NearFieldValue::getZ() const
{
    return m_field(2);
}

bool NearFieldValue::operator==(const NearFieldValue& rhs) const
{
    if ( &rhs == this )
    {
      return true;
    }

    if (rhs.getX() != m_field(0))
    {
        return false;
    }

    if (rhs.getY() != m_field(1))
    {
        return false;
    }

    if (rhs.getZ() != m_field(2))
    {
        return false;
    }

    return true;
}

bool NearFieldValue::operator!=(const NearFieldValue& rhs) const
{
    return not ( this->operator==(rhs) );
}

NearFieldValue& NearFieldValue::operator=(const NearFieldValue& rhs)
{
    if ( &rhs == this )
    {
      return *this;
    }

    m_field(0) = rhs.getX();
    m_field(1) = rhs.getY();
    m_field(2) = rhs.getZ();

    return *this;
}

NearFieldValue NearFieldValue::operator+(const NearFieldValue& rhs) const
{
    NearFieldValue r;

    r.setX(m_field(0) + rhs.getX());
    r.setY(m_field(1) + rhs.getY());
    r.setZ(m_field(2) + rhs.getZ());
    return r;
}

NearFieldValue& NearFieldValue::operator+=(const NearFieldValue& rhs)
{
    *this = *this + rhs;
    return *this;
}

NearFieldValue NearFieldValue::operator-(const NearFieldValue& rhs) const
{
    NearFieldValue r;

    r.setX(m_field(0) - rhs.getX());
    r.setY(m_field(1) - rhs.getY());
    r.setZ(m_field(2) - rhs.getZ());
    return r;
}

NearFieldValue& NearFieldValue::operator-=(const NearFieldValue& rhs)
{
    *this = *this - rhs;
    return *this;
}

// JIF: Need to add tests
bool NearFieldValue::tolerantEqualTo(const NearFieldValue &rhs) const
{
    if (rhs == *this)
    {
        return true;
    }

    double tolerance {1e-15};

    if (std::abs(rhs.getX() - m_field(0)) > tolerance)
    {
        return false;
    }

    if (std::abs(rhs.getY() - m_field(1)) > tolerance)
    {
        return false;
    }

    if (std::abs(rhs.getZ() - m_field(2)) > tolerance)
    {
        return false;
    }

    return true;

}
