#include "MoM.hpp"
#include "assert.h"
#include <iostream>
#include "EdgeContainer.hpp"
#include "Quadrature.hpp"

MoM::MoM(TriangleContainer &tContainer): m_tContainer(tContainer)
{
    //ctor
    m_frequency = -1.0; // Default to an invalid value (should rather be a defined constant)
}

MoM::~MoM()
{
    //dtor
}

void MoM::setFrequency(double freq)
{
    assert(freq > 0.0);
    m_frequency = freq;
}

ComplexMatrix MoM::fillZmatrixTriangle()
{
    return fillZmatrixTriangleInefficient();
}

ComplexMatrix MoM::fillZmatrixTriangleInefficient()
{
    assert(m_frequency > 0.0);

    const std::complex<double> j (0.0, -1.0);
    const double eps0 {8.85418781761e-12};
    const double mu0  {1.25663706143592e-06};
    const double c0 {299792456.2};

    double omega = 2*M_PI*m_frequency;
    double k = (2*M_PI*m_frequency)/c0;

    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();
    std::cout << "Triangle container size: " << m_tContainer.size() << std::endl;
    std::cout << "Edge container size:     " << eContainer.size() << std::endl;

    ComplexMatrix Zmatrix {(unsigned)eContainer.size()};

    // For each edge/basis function m
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        auto triangleList_m = eContainer.at(m).getTriangles();
//        std::cout << triangleList_m.at(0) << "   " << triangleList_m.at(1) << std::endl;

        // Each triangle at the edge (edge always only has two triangles - we don't support junctions)
        for (unsigned tmindex = 0; tmindex < 2 ; ++tmindex)
        {
            Triangle triangle_m = m_tContainer.at(triangleList_m.at(tmindex));
            triangle_m.setOppositeEdge(eContainer.at(m).n1(), eContainer.at(m).n2());
            Node rmc = triangle_m.centre();
//            std::cout << "m( "  << pContainer.find(triangle_m.n1())
//                      << " , " << pContainer.find(triangle_m.n2())
//                      << " , " << pContainer.find(triangle_m.n3())
//                      << " )" << std::endl;

//            std::cout << "( "  << rmc.x()
//                      << " , " << rmc.y()
//                      << " , " << rmc.z()
//                      << " )" << std::endl;

            double sign_m = -(double)tmindex*2.0 + 1.0; // First triangle is positive, second triangle is negative
            assert((sign_m == -1.0) || (sign_m == 1.0));
            Node rho_m_c = (rmc - triangle_m.n1()) * sign_m; // Vector from opposite vertex to centre of triangle (for positive triangle)
//            std::cout << "( "  << rho_m_c.x()
//                      << " , " << rho_m_c.y()
//                      << " , " << rho_m_c.z()
//                      << " )" << std::endl;

            // For each edge/basis function n
            for (unsigned n = 0 ; n < eContainer.size() ; ++n)
            {
                auto triangleList_n = eContainer.at(n).getTriangles();

                // Each triangle at the edge (edge always only has two triangles)
                for (unsigned tnindex = 0 ; tnindex < 2 ; ++tnindex)
//                for (unsigned tnindex = tmindex; tnindex == tmindex ; ++tnindex)
                {
                    Triangle triangle_n = m_tContainer.at(triangleList_n.at(tnindex));
                    triangle_n.setOppositeEdge(eContainer.at(n).n1(), eContainer.at(n).n2());

//                    std::cout << "n( "  << pContainer.find(triangle_n.n1())
//                              << " , " << pContainer.find(triangle_n.n2())
//                              << " , " << pContainer.find(triangle_n.n3())
//                              << " )" << std::endl;


                    double sign_n = -(double)tnindex*2.0 + 1.0; // First triangle is positive, second triangle is negative
                    assert((sign_n == -1.0) || (sign_n == 1.0));

                    // Calculate quadrature points
                    std::vector<Quadrature::WeightedPoint> weightedPoints;
                    weightedPoints = Quadrature::getTriangleSimplexGaussianQuadraturePoints(6);

                    std::complex<double> Phi {0.0, 0.0};
                    std::complex<double> A [3] { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

                    // Loop over qudrature points (to perform integration)
                    for (unsigned qpIndex = 0; qpIndex < weightedPoints.size() ; ++qpIndex)
                    {
                        Quadrature::WeightedPoint wp;
                        wp = weightedPoints.at(qpIndex);
                        Node r = triangle_n.fromSimplex(wp.node);
                        Node rho = sign_n * (r - triangle_n.n1());

//                        std::cout << "wp "  << qpIndex
//                                  << " ( "  << r.x()
//                                  << " , " << r.y()
//                                  << " , " << r.z()
//                                  << " )" << std::endl;

//                        std::cout << "( "  << rho.x()
//                                  << " , " << rho.y()
//                                  << " , " << rho.z()
//                                  << " )" << std::endl;

                        // Calculate Rm
                        double Rm = rmc.distance(r);
//                        std::cout << "( "  << Rm << std::endl;

                        // Calculate scalar potential Phi  (Why Wilco does not have negative here?)
                        Phi += -1.0/(4.0*M_PI*j*omega*eps0) *
                                 ( sign_n * eContainer.at(n).length()/triangle_n.area()) *
                                 (exp(-1.0*j*k*Rm) / Rm) *
                                 (wp.weight*triangle_n.area());
//                        std::cout << "Phi = "  << Phi << std::endl;

                        // Calculate vector potential A
                        std::complex<double> Atmp;
                        Atmp = mu0/(4.0*M_PI) *
                                ( eContainer.at(n).length()/( triangle_n.area() * 2.0) ) *
                                (exp(-1.0*j*k*Rm) / Rm) *
                                (wp.weight*triangle_n.area());
                        A[0] += Atmp * rho.x();
                        A[1] += Atmp * rho.y();
                        A[2] += Atmp * rho.z();
//                        std::cout << "A[0] = "  << A[0] << std::endl;
//                        std::cout << "A[1] = "  << A[1] << std::endl;
//                        std::cout << "A[2] = "  << A[2] << std::endl;
                    }

                    // Add contribution to Z_mn
                    std::complex<double> AdotRho = A[0]*rho_m_c.x() + A[1]*rho_m_c.y() +  A[2]*rho_m_c.z();
//                    std::cout << "AdotRho = "  << AdotRho << std::endl;

                    Zmatrix(m, n) += eContainer.at(m).length() * ( j*omega * (AdotRho/2.0) - sign_n*(Phi) );
                    std::cout << m << ", " << n <<  "  ( " << Zmatrix(m, n) << " " << tmindex << " " << tnindex << "   " << sign_n << " " << sign_m << std::endl;


//					std::complex<double> zero {0.0, 0.0};
//                    if ( Zmatrix(m, n) == zero )
//                    {
//                        std::cout << m << ", " << n << std::endl;
//                        std::cout << "::        " <<
//                                     eContainer.at(m).length() * ( j*omega * (AdotRho/2.0) + (Phi_mn) ) << "   " <<
//                                     std::endl;
//                    }
                }
            }
        }
    }

    return Zmatrix;
}

ComplexMatrix MoM::calculateRHS()
{
    const std::complex<double> j (0.0, -1.0);
//    const double eps0 {8.85418781761e-12};
//    const double mu0  {1.25663706143592e-06};
//    const double c0 {299792456.2};

//    double omega = 2*M_PI*freq;
//    double k = (2*M_PI*freq)/c0;

    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    ComplexMatrix Vvector {(unsigned)eContainer.size(), 1};

    // For each edge/basis function m
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        auto triangleList_m = eContainer.at(m).getTriangles();

        // Each triangle at the edge (edge always only has two triangles - we don't support junctions)
        for (unsigned tmindex = 0; tmindex < 2 ; ++tmindex)
        {
            Triangle triangle_m = m_tContainer.at(triangleList_m.at(tmindex));
            triangle_m.setOppositeEdge(eContainer.at(m).n1(), eContainer.at(m).n2());
            Node rmc = triangle_m.centre();

            double sign_m = -(double)tmindex*2.0 + 1.0; // First triangle is positive, second triangle is negative
            assert((sign_m == -1.0) || (sign_m == 1.0));
            Node rho_m_c = (rmc - triangle_m.n1()) * sign_m; // Vector from opposite vertex to centre of triangle (for positive triangle)

            std::complex<double> Em [3];
            Em[0] = exp(j*rmc.z());
            Em[1] = 0.0;
            Em[2] = 0.0;

            std::complex<double> EdotRho = Em[0]*rho_m_c.x() + Em[1]*rho_m_c.y() +  Em[2]*rho_m_c.z();
            Vvector(m, 0) += eContainer.at(m).length() * ( EdotRho/2.0 );
        }
    }

    return Vvector;
}


// Efficient matrix fill
//    for (TriangleContainer::SizeType observationIndex = 0; observationIndex < tContainer.size() ; ++observationIndex)
//    {
//        Triangle tObserver {tContainer.at(observationIndex)};

//        std::vector<Quadrature::WeightedPoint> weightedObsPointsSimplex;
//        std::vector<Quadrature::WeightedPoint> weightedObsPoints;
//        weightedObsPointsSimplex = Quadrature::getTriangleSimplexGaussianQuadraturePoints(3);
//        weightedObsPoints = weightedObsPointsSimplex;
//        for (unsigned ii = 0; ii < weightedObsPoints.size() ; ++ii)
//        {
//            weightedObsPoints.at(ii).node = tObserver.fromSimplex(weightedObsPointsSimplex.at(ii).node) ;
//        }

//        for (TriangleContainer::SizeType sourceIndex = 0; sourceIndex < tContainer.size() ; ++sourceIndex)
//        {
//            Triangle tSource {tContainer.at(sourceIndex)};

//            std::vector<Quadrature::WeightedPoint> weightedSourcePointsSimplex;
//            std::vector<Quadrature::WeightedPoint> weightedSourcePoints;
//            weightedSourcePointsSimplex = Quadrature::getTriangleSimplexGaussianQuadraturePoints(3);
//            weightedSourcePoints = weightedSourcePointsSimplex;
//            for (unsigned ii = 0; ii < weightedSourcePoints.size() ; ++ii)
//            {
//                weightedSourcePoints.at(ii).node = tSource.fromSimplex(weightedSourcePointsSimplex.at(ii).node) ;
//            }

//            bool requiresSingularIntegration = tContainer.hasCommonNode(observationIndex, sourceIndex);
//            std::cout << requiresSingularIntegration << std::endl;

//        }
//    }

