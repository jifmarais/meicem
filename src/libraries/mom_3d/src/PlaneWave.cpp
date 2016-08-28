#include "PlaneWave.hpp"
#include "math.h"
#include "EMconst.hpp"
#include "assert.h"

PlaneWave::PlaneWave()
{
    //ctor
}

PlaneWave::~PlaneWave()
{
    //dtor
}

void PlaneWave::setAngleOfIncidence(double theta, double phi)
{
    m_incidentTheta = theta/180*EMconst::pi;
    m_incidentPhi   = phi/180*EMconst::pi;
}

void PlaneWave::setPolarisationAngle(double eta)
{
    // Following the convention in FEKO
    // JIF: I'm not really 100% sure this is correct. Need to check this
    // and how FEKO really defines its plane wave
    m_polarisationEta = (180 - eta)/180*EMconst::pi;
}

void PlaneWave::setAmplitude(double amplitude)
{
    m_amplitude = amplitude;
}

void PlaneWave::setFrequency(double frequency)
{
    assert(frequency > 0);
    m_frequency = frequency;
}

void PlaneWave::setFieldPoint(Node p)
{
    m_fieldPoint.set(p.x(), p.y(), p.z());
}

void PlaneWave::setFieldPoint(double x, double y, double z)
{
    m_fieldPoint.set(x, y, z);
}

NearFieldValue PlaneWave::getElectricField()
{
    NearFieldValue val;
    double k = (2*EMconst::pi*m_frequency)/EMconst::c0;

    Vector directionOfIncidence {-sin(m_incidentTheta)*cos(m_incidentPhi) ,
                                 -sin(m_incidentTheta)*sin(m_incidentPhi) ,
                                 -cos(m_incidentTheta)};

    double kDOTr = k*Vector::dot(directionOfIncidence, m_fieldPoint);

    std::complex<double> Etheta = m_amplitude*cos(m_polarisationEta)*exp(EMconst::j*kDOTr);
    std::complex<double> Ephi   = m_amplitude*sin(m_polarisationEta)*exp(EMconst::j*kDOTr);

    std::complex<double> Ex = Etheta*cos(m_incidentTheta)*cos(m_incidentPhi) - Ephi*sin(m_incidentPhi);
    std::complex<double> Ey = Etheta*cos(m_incidentTheta)*sin(m_incidentPhi) + Ephi*cos(m_incidentPhi);
    std::complex<double> Ez = -Etheta*sin(m_incidentTheta);

    val.setX(Ex);
    val.setY(Ey);
    val.setZ(Ez);
    return val;
}

NearFieldValue PlaneWave::getMagneticField()
{
    NearFieldValue Efield = getElectricField();
    NearFieldValue Hfield;
    Vector val;
    val.set(std::abs(Efield.getX()),
            std::abs(Efield.getY()),
            std::abs(Efield.getZ()));
    val = val.norm();
    Vector directionOfIncidence {-sin(m_incidentTheta)*cos(m_incidentPhi) ,
                                 -sin(m_incidentTheta)*sin(m_incidentPhi) ,
                                 -cos(m_incidentTheta)};
    double Zf = 1 / (EMconst::eps0*EMconst::c0);
    Hfield.setX((directionOfIncidence.y()*Efield.getZ() - directionOfIncidence.z()*Efield.getY()) / Zf);
    Hfield.setY((directionOfIncidence.z()*Efield.getX() - directionOfIncidence.x()*Efield.getZ()) / Zf);
    Hfield.setZ((directionOfIncidence.x()*Efield.getY() - directionOfIncidence.y()*Efield.getX()) / Zf);

    return Hfield;
}
