#pragma once
#include <armadillo>
#include <complex>

class NearFieldValue
{
    public:
        typedef arma::cx_vec::fixed<3> cvec;
        NearFieldValue();
        virtual        		~NearFieldValue();
        void			set(std::complex<double> x,
                                    std::complex<double> y,
                                    std::complex<double> z);
        void			setX(std::complex<double> val);
        void			setY(std::complex<double> val);
        void			setZ(std::complex<double> val);
        cvec 			getVec() const;
        void 			setVec(cvec value);
        std::complex<double>	getX() const;
        std::complex<double>	getY() const;
        std::complex<double>	getZ() const;
        NearFieldValue&         operator=(const NearFieldValue& rhs);
        bool                    operator==(const NearFieldValue& rhs) const;
        bool                    operator!=(const NearFieldValue& rhs) const;
        NearFieldValue          operator+(const NearFieldValue& rhs) const;
        NearFieldValue&         operator+=(const NearFieldValue& rhs);
        NearFieldValue          operator-(const NearFieldValue& rhs) const;
        NearFieldValue&         operator-=(const NearFieldValue& rhs);
        bool                    tolerantEqualTo(const NearFieldValue& rhs) const;

    protected:
    private:
        cvec m_field;
};
