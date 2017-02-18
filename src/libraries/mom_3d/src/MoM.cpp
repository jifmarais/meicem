#include "MoM.hpp"
#include "assert.h"
#include <iostream>
#include <iomanip>      // std::setprecision
#include "fstream"
#include "EdgeContainer.hpp"
#include "Quadrature.hpp"
#include "EMconst.hpp"
#include "PlaneWave.hpp"
#include "NearFieldValue.hpp"

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

ComplexMatrix MoM::fillZmatrixTriangleEfficient()
{
    // Here we fill the matrix by looping over the triangles instead of the
    // edges. This is more efficient. Initially, we will be using single point
    // quadrature for the observation triangle and 6 point quadrature for the
    // source triangle (to avoid singularities).

    assert(m_frequency > 0.0);

//    double omega = 2 * EMconst::pi*m_frequency;
//    double k = (2 * EMconst::pi * m_frequency) / EMconst::c0;

    // Calculate quadrature points - we use 6 to avoid sigularity at the centre
    std::vector<Quadrature::WeightedPoint> weightedPointsObservation;
    weightedPointsObservation = Quadrature::getTriangleSimplexGaussianQuadraturePoints(1);
    std::vector<Quadrature::WeightedPoint> weightedPointsSource;
    weightedPointsSource = Quadrature::getTriangleSimplexGaussianQuadraturePoints(6);

    // Container to hold basis functions (non-boundary edges)
    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    ComplexMatrix Zmatrix {(unsigned)eContainer.size()};

    // For each observation triangle
    for (unsigned toIndex = 0 ; toIndex < m_tContainer.size() ; ++toIndex)
    {
        // For each source triangle
        for (unsigned toIndex = 0 ; toIndex < m_tContainer.size() ; ++toIndex)
        {
            // For each non-boundary edge of the observation triangle (basis function)
            //getEdgeIndecesOnTriangle
            for (unsigned toEdgeIndex = 0 ; toEdgeIndex < m_tContainer.size() ; ++toEdgeIndex)
            {

            }

        }
    }

    return Zmatrix;
}

std::complex<double> MoM::G0(const double R, const double k) const
{
    return exp(-1.0 * EMconst::j * k * R) / (4.0 * EMconst::pi * R);
}

double MoM::RWGBasisFunction(const Triangle T) const
{
    double ln = T.n2().distance(T.n3());
    double A  = T.area();
    //return (ln/(2*A))*simplexCoordinate;
    return (ln/(2.0*A));
}

double MoM::divRWGBasisFunction(const Triangle T, const double sign) const
{
    double ln = T.n2().distance(T.n3());
    double A  = T.area();
    return sign*ln/A;
}

ComplexMatrix MoM::fillZmatrixTriangleInefficient1()
{
    // Here we fill the matrix and simply use 1 point integration at the
    // observation triangle and 6 at the source triangle to avoid singularities.
    // Not a good idea, but it gets things going.

    assert(m_frequency > 0.0);

    double omega = 2 * EMconst::pi*m_frequency;
    double k = (2 * EMconst::pi * m_frequency) / EMconst::c0;

    // Calculate quadrature points - we use 6 to avoid sigularity at the centre
    std::vector<Quadrature::WeightedPoint> weightedPoints;
    weightedPoints = Quadrature::getTriangleSimplexGaussianQuadraturePoints(6);

    // Container to hold basis functions (non-boundary edges)
    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    ComplexMatrix Zmatrix {(unsigned)eContainer.size()};

    // For each edge/basis function m
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        // Triangles bounding the current edge (index m)
        auto edge_m = eContainer.at(m);
        auto boundingTrianglesIndexList_m = edge_m.getTriangles();

        // Each triangle at the edge (edge always only has two triangles - we don't support junctions yet)
        assert(boundingTrianglesIndexList_m.size() == 2);
        for (unsigned triangle_mIndex = 0; triangle_mIndex < boundingTrianglesIndexList_m.size() ; ++triangle_mIndex)
        {
            Triangle triangle_m = m_tContainer.at(boundingTrianglesIndexList_m.at(triangle_mIndex));

            //set triangle n1 to be the vertex opposite to the edge
            triangle_m.setOppositeEdge(edge_m.n1(), edge_m.n2());

            Node r_mCentre = triangle_m.centre();

            // Triangle current is from triangle with smaller index to larger index
            // First triangle is positive, second triangle is negative (they are sorted)
            double sign_m = -(double)triangle_mIndex*2.0 + 1.0;
            assert((sign_m == -1.0) || (sign_m == 1.0));

            // Vector from opposite vertex to centre of triangle (for positive triangle)
            Node rho_mCentre = (r_mCentre - triangle_m.n1()) * sign_m;

            // For each edge/basis function n
            for (unsigned n = 0 ; n < eContainer.size() ; ++n)
            {
                // Triangles bounding the current edge (index n)
                auto edge_n = eContainer.at(n);
                auto boundingTrianglesIndexList_n = edge_n.getTriangles();

                // Each triangle at the edge (edge always only has two triangles since we don't support junctions yet)
                assert(boundingTrianglesIndexList_n.size() == 2);
                for (unsigned triangle_nIndex = 0 ; triangle_nIndex < boundingTrianglesIndexList_n.size() ; ++triangle_nIndex)
                {
                    Triangle triangle_n = m_tContainer.at(boundingTrianglesIndexList_n.at(triangle_nIndex));

                    //set triangle n1 to be the vertex opposite to the edge
                    triangle_n.setOppositeEdge(edge_n.n1(), edge_n.n2());

                    // First triangle is positive, second triangle is negative
                    double sign_n = -(double)triangle_nIndex*2.0 + 1.0;
                    assert((sign_n == -1.0) || (sign_n == 1.0));

                    std::complex<double> expValue {0.0, 0.0};
                    std::complex<double> Phi {0.0, 0.0};
                    std::complex<double> A [3] { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

                    // Loop over qudrature points (to perform integration)
                    for (unsigned qpIndex = 0; qpIndex < weightedPoints.size() ; ++qpIndex)
                    {
                        Quadrature::WeightedPoint wp;
                        wp = weightedPoints.at(qpIndex);

                        // Vector to point from the origin
                        Node r_n_qp = triangle_n.fromSimplex(wp.node);
                        Node rho_n_qp = sign_n * (r_n_qp - triangle_n.n1());

                        // Calculate Rm
                        double Rm = r_mCentre.distance(r_n_qp);

                        // Calculate scalar potential Phi
                        expValue =  G0(Rm, k) * wp.weight;
                        Phi += expValue;

                        // Calculate vector potential A
                        A[0] += expValue * rho_n_qp.x();
                        A[1] += expValue * rho_n_qp.y();
                        A[2] += expValue * rho_n_qp.z();
                   }
                    Phi  *= -1.0 / (EMconst::j * omega * EMconst::eps0) * ( sign_n * edge_n.length() );

                    std::complex<double> Atmp;
                    Atmp  = EMconst::mu0 * ( edge_n.length() / 2.0 ) ;

                    // Add contribution to Z_mn
                    std::complex<double> AdotRho = Atmp * ( A[0]*rho_mCentre.x() + A[1]*rho_mCentre.y() +  A[2]*rho_mCentre.z() );

                    Zmatrix(m, n) += edge_m.length() * ( EMconst::j * omega * (AdotRho/2.0) - sign_m*(Phi) );
                }
            }
        }
    }

    return Zmatrix;
}

ComplexMatrix MoM::fillZmatrixTriangleInefficient()
{
    // Here we fill the matrix and simply use 1 point integration at the
    // observation triangle and 6 at the source triangle to avoid singularities.
    // Not a good idea, but it gets things going.

    assert(m_frequency > 0.0);

    //double omega = 2 * EMconst::pi*m_frequency;
    double k = (2 * EMconst::pi * m_frequency) / EMconst::c0;
    double accurateIntegrationDistance = (EMconst::c0 / m_frequency) / 15.0 ;
    double aiCount = 0;

    // Calculate quadrature points - we use 6 to avoid sigularity at the centre
    std::vector<Quadrature::WeightedPoint> weightedPointsSource = Quadrature::getTriangleSimplexGaussianQuadraturePoints(6);
    std::vector<Quadrature::WeightedPoint> weightedPointsTest   = Quadrature::getTriangleSimplexGaussianQuadraturePoints(3);

    // Container to hold basis functions (non-boundary edges)
    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    ComplexMatrix Zmatrix {(unsigned)eContainer.size()};

    // For each edge/basis function m
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        // Triangles bounding the current edge (index m)
        auto edge_m = eContainer.at(m);
        auto boundingTrianglesIndexList_m = edge_m.getTriangles();

        // Each triangle at the edge (edge always only has two triangles - we don't support junctions yet)
        assert(boundingTrianglesIndexList_m.size() == 2);
        for (unsigned triangle_mIndex = 0; triangle_mIndex < boundingTrianglesIndexList_m.size() ; ++triangle_mIndex)
        {
            Triangle triangle_m = m_tContainer.at(boundingTrianglesIndexList_m.at(triangle_mIndex));

            //set triangle n1 to be the vertex opposite to the edge
            triangle_m.setOppositeEdge(edge_m.n1(), edge_m.n2());

            // Triangle current is from triangle with smaller index to larger index
            // First triangle is positive, second triangle is negative (they are sorted)
            double sign_m = -(double)triangle_mIndex*2.0 + 1.0;
            assert((sign_m == -1.0) || (sign_m == 1.0));

            // Loop over qudrature points (to perform outer integration)
            for (unsigned qpmIndex = 0; qpmIndex < weightedPointsTest.size() ; ++qpmIndex)
            {
                Quadrature::WeightedPoint wpm = weightedPointsTest.at(qpmIndex);

                // Vector to point from the origin
                Node r_m_qp = triangle_m.fromSimplex(wpm.node);
                Node rho_m_qp = sign_m * (r_m_qp - triangle_m.n1());

                // For each edge/basis function n
                for (unsigned n = 0 ; n < eContainer.size() ; ++n)
                {
                    // Triangles bounding the current edge (index n)
                    auto edge_n = eContainer.at(n);
                    auto boundingTrianglesIndexList_n = edge_n.getTriangles();

                    // Each triangle at the edge (edge always only has two triangles since we don't support junctions yet)
                    assert(boundingTrianglesIndexList_n.size() == 2);
                    for (unsigned triangle_nIndex = 0 ; triangle_nIndex < boundingTrianglesIndexList_n.size() ; ++triangle_nIndex)
                    {
                        Triangle triangle_n = m_tContainer.at(boundingTrianglesIndexList_n.at(triangle_nIndex));

                        //set triangle n1 to be the vertex opposite to the edge
                        triangle_n.setOppositeEdge(edge_n.n1(), edge_n.n2());

                        // First triangle is positive, second triangle is negative
                        double sign_n = -(double)triangle_nIndex*2.0 + 1.0;
                        assert((sign_n == -1.0) || (sign_n == 1.0));


                        //if (r_m_qp.distance(triangle_n.centre()) <  accurateIntegrationDistance)
                        if (triangle_m.centre().distance(triangle_n.centre()) <  accurateIntegrationDistance)
                        {
                            aiCount += 1;
                        }

                        std::complex<double> A [3] { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
                        std::complex<double> B {0.0, 0.0};

                        // Loop over qudrature points (to perform integration)
                        for (unsigned qpnIndex = 0; qpnIndex < weightedPointsSource.size() ; ++qpnIndex)
                        {
                            Quadrature::WeightedPoint wpn = weightedPointsSource.at(qpnIndex);

                            // Vector to point from the origin
                            Node r_n_qp = triangle_n.fromSimplex(wpn.node);
                            Node rho_n_qp = sign_n * (r_n_qp - triangle_n.n1());

                            double R = r_m_qp.distance(r_n_qp);
                            std::complex<double> tmp = G0(R, k) * wpn.weight;

                            // Calculate A
                            A[0] += rho_n_qp.x() * tmp;
                            A[1] += rho_n_qp.y() * tmp;
                            A[2] += rho_n_qp.z() * tmp;

                            // Calculate B
                            B += tmp;
                        }
                        B *= divRWGBasisFunction(triangle_n, sign_n);
                        A[0] *= RWGBasisFunction(triangle_n);
                        A[1] *= RWGBasisFunction(triangle_n);
                        A[2] *= RWGBasisFunction(triangle_n);
                        std::complex<double> Adotfm = RWGBasisFunction(triangle_m) * ( A[0]*rho_m_qp.x() + A[1]*rho_m_qp.y() +  A[2]*rho_m_qp.z() ) * wpm.weight;
                        std::complex<double> tmp = divRWGBasisFunction(triangle_m, sign_m) * B * wpm.weight;

                        Zmatrix(m, n) += EMconst::j * EMconst::Z0 * (k*Adotfm - tmp/k) * triangle_m.area()*triangle_n.area();
                    }
                }
            }
        }
    }
    std::cout << "accurate: " << aiCount << std::endl;

    return Zmatrix;
}

ComplexMatrix MoM::calculateRHS()
{
    PlaneWave pw;
    pw.setFrequency(m_frequency);

    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    ComplexMatrix Vvector {(unsigned)eContainer.size(), 1};
    std::vector<Quadrature::WeightedPoint> weightedPointsTest   = Quadrature::getTriangleSimplexGaussianQuadraturePoints(1);

    // For each edge/basis function m
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        auto triangleList_m = eContainer.at(m).getTriangles();

        // Each triangle at the edge (edge always only has two triangles - we don't support junctions)
        for (unsigned tmindex = 0; tmindex < 2 ; ++tmindex)
        {
            Triangle triangle_m = m_tContainer.at(triangleList_m.at(tmindex));
            triangle_m.setOppositeEdge(eContainer.at(m).n1(), eContainer.at(m).n2()); // Set the remaining vertex as the first node of the triangle
            //Node rmc = triangle_m.centre();

            double sign_m = -(double)tmindex*2.0 + 1.0; // First triangle is positive, second triangle is negative
            assert((sign_m == -1.0) || (sign_m == 1.0));

            for (unsigned qpmIndex = 0; qpmIndex < weightedPointsTest.size() ; ++qpmIndex)
            {
                Quadrature::WeightedPoint wpm = weightedPointsTest.at(qpmIndex);

                // Vector to point from the origin
                Node r_m_qp = triangle_m.fromSimplex(wpm.node);
                Node rho_m_qp = sign_m * (r_m_qp - triangle_m.n1());

                pw.setFieldPoint(r_m_qp);
                NearFieldValue Em = pw.getElectricField();
                std::complex<double> Edotfm = RWGBasisFunction(triangle_m) * (Em.getX()*rho_m_qp.x() + Em.getY()*rho_m_qp.y() +  Em.getZ()*rho_m_qp.z()) * wpm.weight;
                Vvector(m, 0) += triangle_m.area()*Edotfm;
           }
        }
    }

    return Vvector;
}

void MoM::writeCurrentsToOS(std::string fname, ComplexMatrix solutionMatrix) const
{
    // JIF: This is the wrong place for the method, but I just want to get it working. Should be refactored.

    assert(fname != "");

    std::string ext {"os"};
    std::ofstream file (fname + '.' + ext);

    if (file.is_open())
    {
        file << "##File Type: Currents" << std::endl;
        file << "##File Format: 4" << std::endl;
        file << "** File exported by MEICEM" << std::endl;
        file << std::endl;
//        file << "#Configuration Name: StandardConfiguration1" << std::endl;
//        file << "#Request Name: NearField11" << std::endl;
        file << "#Frequency:   1.00000000E+008" << std::endl;
        file << "#No. of Electric Current Triangle Samples: " << m_tContainer.size() << std::endl;
        file << "#No. of Header Lines: 1" << std::endl;
        file << "#";
        file << std::setw(10);
        file << "\"Num\"";
        file << std::setw(16);
        file << "\"X\"";
        file << std::setw(19);
        file << "\"Y\"";
        file << std::setw(19);
        file << "\"Z\"";
        file << std::setw(19);
        file << "\"Re(Jx)\"";
        file << std::setw(19);
        file << "\"Im(Jx)\"";
        file << std::setw(19);
        file << "\"Re(Jy)\"";
        file << std::setw(19);
        file << "\"Im(Jy)\"";
        file << std::setw(19);
        file << "\"Re(Jz)\"";
        file << std::setw(19);
        file << "\"Im(Jz)\"";
        file << std::setw(19);
        file << "\"Abs(Jcorn1)\"";
        file << std::setw(19);
        file << "\"Abs(Jcorn2)\"";
        file << std::setw(19);
        file << "\"Abs(Jcorn3)\"";
        file << std::setw(19);
        file << "\"Re(Jx_c1)\"";
        file << std::setw(19);
        file << "\"Im(Jx_c1)\"";
        file << std::setw(19);
        file << "\"Re(Jy_c1)\"";
        file << std::setw(19);
        file << "\"Im(Jy_c1)\"";
        file << std::setw(19);
        file << "\"Re(Jz_c1)\"";
        file << std::setw(19);
        file << "\"Im(Jz_c1)\"";
        file << std::setw(19);
        file << "\"Re(Jx_c2)\"";
        file << std::setw(19);
        file << "\"Im(Jx_c2)\"";
        file << std::setw(19);
        file << "\"Re(Jy_c2)\"";
        file << std::setw(19);
        file << "\"Im(Jy_c2)\"";
        file << std::setw(19);
        file << "\"Re(Jz_c2)\"";
        file << std::setw(19);
        file << "\"Im(Jz_c2)\"";
        file << std::setw(19);
        file << "\"Re(Jx_c3)\"";
        file << std::setw(19);
        file << "\"Im(Jx_c3)\"";
        file << std::setw(19);
        file << "\"Re(Jy_c3)\"";
        file << std::setw(19);
        file << "\"Im(Jy_c3)\"";
        file << std::setw(19);
        file << "\"Re(Jz_c3)\"";
        file << std::setw(19);
        file << "\"Im(Jz_c3)\"";
        file << std::endl;

        EdgeContainer eContainer(m_tContainer);
        eContainer.buildNonboundaryEdgeList();

        for (unsigned tIndex=0; tIndex < m_tContainer.size(); ++tIndex)
        {
//            file << std::fixed ;
            file << std::setprecision(8) ;
            file << std::scientific ;
            file << std::setw(10);
            file << tIndex + 1 << " ";
            file << std::setw(21);
            file << m_tContainer.at(tIndex).centre().x();
            file << std::setw(19);
            file << m_tContainer.at(tIndex).centre().y();
            file << std::setw(19);
            file << m_tContainer.at(tIndex).centre().z();

            // Use near field values to store currents (complex vectors)
            std::vector<NearFieldValue> nodeCurrents;
            for (auto nIndex=0; nIndex < 3 ; ++nIndex)
            {
                NearFieldValue val;
                nodeCurrents.push_back(val);
            }

            // Find the non-boundary edges for a given triangle
            std::vector<EdgeContainer::SizeType> edgeIndices;
            edgeIndices = eContainer.getEdgeIndecesOnTriangle(tIndex);

//            std::cout << "Triangle " << tIndex+1 << ": " << std::endl;

            for (auto nIndex=0; nIndex < 3 ; ++nIndex)
            {
                Node n;
                n = m_tContainer.at(tIndex)[nIndex];
//                std::cout << "  n" << nIndex << ": " << n.x() << " , " << n.y() << " , " << n.z() << std::endl;
                Node tNormal;
                tNormal = m_tContainer.at(tIndex).normal();



                for (EdgeContainer::SizeType eIndex=0; eIndex < edgeIndices.size() ; ++eIndex)
                {
                    if (eContainer.at(edgeIndices.at(eIndex)).n1() == n ||
                        eContainer.at(edgeIndices.at(eIndex)).n2() == n )
                    {
                        // Calculate edge direction (current normal to the edge)
                        Edge e;
                        e = eContainer.at(edgeIndices.at(eIndex));
                        e.correctOrientation(m_tContainer.at(tIndex));

                        Node otherEdge;
                        for (auto oIndex = 0; oIndex < 3; ++oIndex)
                        {
                            if (m_tContainer.at(tIndex)[oIndex] != e.n1() &&
                                m_tContainer.at(tIndex)[oIndex] != e.n2() )
                            {
                                otherEdge = m_tContainer.at(tIndex)[oIndex] - n;
                                break;
                            }
                        }
                        otherEdge = otherEdge.norm();

                        Node eNormal;
                        eNormal = e.normal();
//                        std::cout << "        edge" << edgeIndices.at(eIndex) << ": " << eNormal.x() << " , " << eNormal.y() << " , " << eNormal.z() << std::endl;
//                        std::cout << "      O edge" << edgeIndices.at(eIndex) << ": " << otherEdge.x() << " , " << otherEdge.y() << " , " << otherEdge.z() << std::endl;

                        Node current;
                        current = Node::cross(tNormal, eNormal);

                        if (tIndex == e.getTriangles().at(0))
                        {
                            // Current flowing out of this triangle
                            current *= -1.0;
                        }
                        current = otherEdge*Node::dot(current, otherEdge);

//                        std::cout << "              current: " << nIndex << " | " << current.x() << " , " << current.y() << " , " << current.z() << std::endl;


                        // Scale by coefficient and add not nodal current value
//                        current *= solutionMatrix(edgeIndices.at(eIndex), 1);
                        NearFieldValue val;
                        val = nodeCurrents[nIndex];
//                        std::cout << "val.getY() 1: " << val.getY() << std::endl;
                        val.setX(val.getX() + current.x()*solutionMatrix(edgeIndices.at(eIndex), 0));
                        val.setY(val.getY() + current.y()*solutionMatrix(edgeIndices.at(eIndex), 0));
                        val.setZ(val.getZ() + current.z()*solutionMatrix(edgeIndices.at(eIndex), 0));
                        nodeCurrents[nIndex] = val;
//                        std::cout << "              val    : " << nIndex << " | " << val.getX() << " , " << val.getY() << " , " << val.getZ() << std::endl;

//                        std::cout << "val.getY() 2: " << val.getY() << std::endl;

//                        file << edgeIndices.at(eIndex) << "  ";
                    }
                }
            }

            // Write out current at the node
            file << std::setw(19);
            file << (nodeCurrents.at(0).getX().real() + nodeCurrents.at(1).getX().real() + nodeCurrents.at(2).getX().real()) / 3;
            file << std::setw(19);
            file << (nodeCurrents.at(0).getX().imag() + nodeCurrents.at(1).getX().imag() + nodeCurrents.at(2).getX().imag()) / 3;
            file << std::setw(19);
            file << (nodeCurrents.at(0).getY().real() + nodeCurrents.at(1).getY().real() + nodeCurrents.at(2).getY().real()) / 3;
            file << std::setw(19);
            file << (nodeCurrents.at(0).getY().imag() + nodeCurrents.at(1).getY().imag() + nodeCurrents.at(2).getY().imag()) / 3;
            file << std::setw(19);
            file << (nodeCurrents.at(0).getZ().real() + nodeCurrents.at(1).getZ().real() + nodeCurrents.at(2).getZ().real()) / 3;
            file << std::setw(19);
            file << (nodeCurrents.at(0).getZ().imag() + nodeCurrents.at(1).getZ().imag() + nodeCurrents.at(2).getZ().imag()) / 3;

            file << std::setw(19);
            file << 0.0;
            file << std::setw(19);
            file << 0.0;
            file << std::setw(19);
            file << 0.0;

            for (auto nIndex=0; nIndex < 3 ; ++nIndex)
            {

                file << std::setw(19);
                file << nodeCurrents.at(nIndex).getX().real();
                file << std::setw(19);
                file << nodeCurrents.at(nIndex).getX().imag();
                file << std::setw(19);
                file << nodeCurrents.at(nIndex).getY().real();
                file << std::setw(19);
                file << nodeCurrents.at(nIndex).getY().imag();
                file << std::setw(19);
                file << nodeCurrents.at(nIndex).getZ().real();
                file << std::setw(19);
                file << nodeCurrents.at(nIndex).getZ().imag();

            }


//            file << edgeIndices.size() << "  ";


//            file << m_tContainer.at(tIndex).n1() << "  ";

            file << std::endl;
        }

//        for (unsigned zIndex=0; zIndex < m_numPointsZ ; ++zIndex)
//        {
//            for (unsigned yIndex=0; yIndex < m_numPointsY ; ++yIndex)
//            {
//                for (unsigned xIndex=0; xIndex < m_numPointsX ; ++xIndex)
//                {
//                    file << std::fixed << std::setprecision(9) ;
//                    file << getPointAt(xIndex, yIndex, zIndex).x() << "  " ;
//                    file << getPointAt(xIndex, yIndex, zIndex).y() << "  " ;
//                    file << getPointAt(xIndex, yIndex, zIndex).z() << "  " ;
//                    file << getValueAt(xIndex, yIndex, zIndex).getX().real() << "  " ;
//                    file << getValueAt(xIndex, yIndex, zIndex).getX().imag() << "  " ;
//                    file << getValueAt(xIndex, yIndex, zIndex).getY().real() << "  " ;
//                    file << getValueAt(xIndex, yIndex, zIndex).getY().imag() << "  " ;
//                    file << getValueAt(xIndex, yIndex, zIndex).getZ().real() << "  " ;
//                    file << getValueAt(xIndex, yIndex, zIndex).getZ().imag() ;
//                    file << std::endl;
//                }
//            }
//        }

        file << std::endl;
        file.close();
    }
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

