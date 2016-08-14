#pragma once
#include "NearFieldValue.hpp"
#include "Vector.hpp"
#include <complex>

class PlaneWave
{
    public:
        PlaneWave();
        virtual        		~PlaneWave();
        void			setAngleOfIncidence(double theta, double phi);
        void			setPolarisationAngle(double eta);
        void			setAmplitude(double amplitude);
        void 			setFrequency(double frequency);

        void			setFieldPoint(double x, double y, double z);
        NearFieldValue		getElectricField();
        NearFieldValue		getMagneticField();

protected:
    private:
        double 			m_incidentTheta   {0.0};
        double 			m_incidentPhi     {0.0};
        double 			m_polarisationEta {0.0};
        double			m_amplitude       {1.0};
        double			m_frequency       {1.0};
        Vector			m_fieldPoint;
};
