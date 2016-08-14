#include "NearFieldValue.hpp"
#include "assert.h"

NearFieldValue::NearFieldValue()
{
    //ctor
    m_Xcomponent = {0.0, 0.0};
    m_Ycomponent = {0.0, 0.0};
    m_Zcomponent = {0.0, 0.0};
}

NearFieldValue::~NearFieldValue()
{
    //dtor
}

void NearFieldValue::set(std::complex<double> x, std::complex<double> y, std::complex<double> z)
{
    m_Xcomponent = x;
    m_Ycomponent = y;
    m_Zcomponent = z;
}

void NearFieldValue::setX(std::complex<double> val)
{
    m_Xcomponent = val;
}

void NearFieldValue::setY(std::complex<double> val)
{
    m_Ycomponent = val;
}

void NearFieldValue::setZ(std::complex<double> val)
{
    m_Zcomponent = val;
}

std::complex<double> NearFieldValue::getX() const
{
    return m_Xcomponent;
}

std::complex<double> NearFieldValue::getY() const
{
    return m_Ycomponent;
}

std::complex<double> NearFieldValue::getZ() const
{
    return m_Zcomponent;
}

bool NearFieldValue::operator==(const NearFieldValue& rhs) const
{
    if ( &rhs == this )
    {
      return true;
    }

    if (rhs.getX() != m_Xcomponent)
    {
        return false;
    }

    if (rhs.getY() != m_Ycomponent)
    {
        return false;
    }

    if (rhs.getZ() != m_Zcomponent)
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

    m_Xcomponent = rhs.getX();
    m_Ycomponent = rhs.getY();
    m_Zcomponent = rhs.getZ();

    return *this;
}

NearFieldValue NearFieldValue::operator+(const NearFieldValue& rhs) const
{
    NearFieldValue r;

    r.setX(m_Xcomponent + rhs.getX());
    r.setY(m_Ycomponent + rhs.getY());
    r.setZ(m_Zcomponent + rhs.getZ());
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

    r.setX(m_Xcomponent - rhs.getX());
    r.setY(m_Ycomponent - rhs.getY());
    r.setZ(m_Zcomponent - rhs.getZ());
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

    if (std::abs(rhs.getX() - m_Xcomponent) > tolerance)
    {
        return false;
    }

    if (std::abs(rhs.getY() - m_Ycomponent) > tolerance)
    {
        return false;
    }

    if (std::abs(rhs.getZ() - m_Zcomponent) > tolerance)
    {
        return false;
    }

    return true;

}
