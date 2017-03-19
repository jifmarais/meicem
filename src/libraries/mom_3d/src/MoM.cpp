#include "MoM.hpp"
#include "assert.h"
#include <iostream>
#include <iomanip>      // std::setprecision
#include "fstream"
#include "EdgeContainer.hpp"
#include "Quadrature.hpp"
#include "EMconst.hpp"
#include "NearFieldValue.hpp"

MoM::MoM(TriangleContainer &tContainer): m_tContainer(tContainer)
{
    //ctor
    m_frequency = -1.0; // Default to an invalid value (should rather be a defined constant)
    m_numberOfSourceIntegrationPoints = 3;
    m_numberOfTestIntegrationPoints = 3;
}

MoM::~MoM()
{
    //dtor
}

void MoM::setFrequency(double freq)
{
    assert(freq > 0.0);
    m_frequency = freq;
    m_omega = 2 * EMconst::pi * m_frequency;
    m_k = m_omega / EMconst::c0;
}

void MoM::setNumberOfSourceIntegrationPoints(unsigned n)
{
    assert(n >= 1);
    m_numberOfSourceIntegrationPoints = n;
}

void MoM::setNumberOfTestIntegrationPoints(unsigned n)
{
    assert(n >= 1);
    m_numberOfTestIntegrationPoints = n;
}

arma::cx_mat MoM::fillZmatrixTriangle()
{
//    return fillZmatrixTriangleInefficient();
    return fillZmatrixTriangleEfficient();
}

double MoM::RWGBasisFunctionSign(unsigned baseIndex, unsigned ohterTriangleIndex) const
{
    assert(baseIndex != ohterTriangleIndex);
    int diff  = ohterTriangleIndex - baseIndex;
    return diff/std::fabs(diff);
}

arma::cx_mat MoM::fillZmatrixTriangleEfficient()
{
    // Here we fill the matrix by looping over the triangles instead of the
    // edges. This is more efficient.

    assert(m_frequency > 0.0);

    double accurateIntegrationDistance = (EMconst::c0 / m_frequency) / 5.0 ;

    Quadrature quadrature;

    // Container to hold basis functions (non-boundary edges)
    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    arma::cx_mat Zmatrix ((unsigned)eContainer.size(), (unsigned)eContainer.size());
    Zmatrix.fill(std::complex<double> {0.0, 0.0});

    arma::cx_vec A (3);
    std::complex<double> B {0.0, 0.0};
    std::complex<double> G0weight;
    Quadrature::WeightedPoint wpTest;
    Quadrature::WeightedPoint wpSrc;
    Node r_Src_qp;
    Node::vec rho_Src_qp;
    Node::vec rho_Test_qp;

    // For each test triangle
    for (unsigned tTestIndex = 0 ; tTestIndex < m_tContainer.size() ; ++tTestIndex)
    {
        Triangle testTriangle = m_tContainer.at(tTestIndex);

        // Find the basis functions (edges in this case) that involve this triangle
        std::vector<EdgeContainer::SizeType> testTriangleEdgeList = eContainer.getEdgeIndecesOnTriangle(tTestIndex);

        std::vector<Quadrature::WeightedPoint> weightedPointsTest  = quadrature.getTriangleGaussianQuadraturePoints(testTriangle, m_numberOfTestIntegrationPoints);

        for (unsigned tSrcIndex = 0 ; tSrcIndex < m_tContainer.size() ; ++tSrcIndex)
        {
            Triangle srcTriangle = m_tContainer.at(tSrcIndex);
            std::vector<EdgeContainer::SizeType> sourceTriangleEdgeList = eContainer.getEdgeIndecesOnTriangle(tSrcIndex);
            std::vector<Quadrature::WeightedPoint> defaultWeightedPointsSrc = quadrature.getTriangleGaussianQuadraturePoints(srcTriangle, m_numberOfSourceIntegrationPoints);

            // Loop over qudrature points (to perform outer integration)
            for (unsigned qpTestIndex = 0; qpTestIndex < weightedPointsTest.size() ; ++qpTestIndex)
            {
                wpTest = weightedPointsTest.at(qpTestIndex);

                // Vector to point from the origin
                Node r_Test_qp = wpTest.node;


                std::vector<Quadrature::WeightedPoint> weightedPointsSrc;
                if (srcTriangle.centre().distance(r_Test_qp) <  accurateIntegrationDistance)
                {
                    weightedPointsSrc = quadrature.RAR1S(srcTriangle, r_Test_qp, m_numberOfSourceIntegrationPoints*2);
                }
                else
                {
                    weightedPointsSrc = defaultWeightedPointsSrc;
                }

                //Add matrix contributions
                // For each basis function involved in the test triangle
                for (unsigned testTriangleEdgeIndex = 0; testTriangleEdgeIndex < testTriangleEdgeList.size() ; ++testTriangleEdgeIndex)
                {
                    auto testIndex = testTriangleEdgeList.at(testTriangleEdgeIndex);
                    auto testTriangleEdge = eContainer.at(testIndex);
                    auto trianglesBoundingTestEdge = testTriangleEdge.getTriangles();
                    Node nodeOppositeToTestEdge = testTriangle.getOppositeNode(testTriangleEdge.n1(), testTriangleEdge.n2());
                    rho_Test_qp = r_Test_qp.getVec() - nodeOppositeToTestEdge.getVec(); // Sign handled later
                    double RWGTest = RWGBasisFunction(testTriangle, testTriangleEdge);

                    for (unsigned testTriangleIndex = 0; testTriangleIndex < trianglesBoundingTestEdge.size() ; ++testTriangleIndex)
                    {
                        if (trianglesBoundingTestEdge.at(testTriangleIndex) != tTestIndex)
                        {
                            double signTest = RWGBasisFunctionSign(tTestIndex, trianglesBoundingTestEdge.at(testTriangleIndex));
                            double divRWGTest = divRWGBasisFunction(testTriangle, testTriangleEdge, signTest);

                            // For each basis function involved in the source triangle
                            for (unsigned sourceTriangleEdgeIndex = 0; sourceTriangleEdgeIndex < sourceTriangleEdgeList.size() ; ++sourceTriangleEdgeIndex)
                            {
                                auto sourceIndex = sourceTriangleEdgeList.at(sourceTriangleEdgeIndex);
                                auto srcEdge = eContainer.at(sourceIndex);
                                auto trianglesBoundingSourceEdge = srcEdge.getTriangles();
                                Node nodeOppositeToEdgeSrc = srcTriangle.getOppositeNode(srcEdge.n1(), srcEdge.n2());
                                double RWGSrc = RWGBasisFunction(srcTriangle, srcEdge);
                                double divRWGSrcNoSign = divRWGBasisFunction(srcTriangle, srcEdge, 1.0);

                                for (unsigned srcTriangleIndex = 0; srcTriangleIndex < trianglesBoundingSourceEdge.size() ; ++srcTriangleIndex)
                                {
                                    if (trianglesBoundingSourceEdge.at(srcTriangleIndex) != tSrcIndex)
                                    {
                                        A.fill({0, 0});
                                        B = {0.0, 0.0};

                                        for (unsigned qpSrcIndex = 0; qpSrcIndex < weightedPointsSrc.size() ; ++qpSrcIndex)
                                        {
                                            wpSrc = weightedPointsSrc.at(qpSrcIndex);

                                            // Vector to point from the origin
                                            r_Src_qp = wpSrc.node;
                                            rho_Src_qp = r_Src_qp.getVec() - nodeOppositeToEdgeSrc.getVec(); // Sign handled later

                                            G0weight = G0(r_Test_qp.distance(r_Src_qp)) * wpSrc.weight;

                                            // CRC: Should use real names for quantities (vector potential and scalar potential)
                                            // Calculate A
                                            A += rho_Src_qp * G0weight;

                                            // Calculate B
                                            B += G0weight;
                                        }

                                        double signSrc = RWGBasisFunctionSign(tSrcIndex, trianglesBoundingSourceEdge.at(srcTriangleIndex));

                                        std::complex<double> Adotfm = RWGTest * RWGSrc * signSrc * ( dot(A, rho_Test_qp*signTest) ) * wpTest.weight;
                                        std::complex<double> tmp = divRWGTest * B * divRWGSrcNoSign*signSrc * wpTest.weight;

                                        Zmatrix(testIndex, sourceIndex) += (m_k*Adotfm - tmp/m_k);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Zmatrix *= EMconst::j * EMconst::Z0;

    return Zmatrix;
}

std::complex<double> MoM::G0(const double R) const
{
    return exp(EMconst::j * m_k * R) / (4.0 * EMconst::pi * R);
}

double MoM::RWGBasisFunction(const Triangle T, const Edge E) const
{
    return (E.length()/(2.0*T.area()));
}

double MoM::divRWGBasisFunction(const Triangle T, const Edge E, const double sign) const
{
    return sign*E.length()/T.area();
}

arma::cx_mat MoM::fillZmatrixTriangleInefficient()
{
    // Here we fill the matrix and simply use 1 point integration at the
    // test triangle and 6 at the source triangle to avoid singularities.
    // Not a good idea, but it gets things going.

    assert(m_frequency > 0.0);

    double accurateIntegrationDistance = (EMconst::c0 / m_frequency) / 5.0 ;

    Quadrature quadrature;

    // Container to hold basis functions (non-boundary edges)
    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    arma::cx_mat Zmatrix ((unsigned)eContainer.size(), (unsigned)eContainer.size());
    Zmatrix.fill(std::complex<double> {0.0, 0.0});

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
            auto tmIndex = boundingTrianglesIndexList_m.at(triangle_mIndex);
            Triangle triangle_m = m_tContainer.at(tmIndex);

            //set triangle n1 to be the vertex opposite to the edge
            Node nodeOppositeToEdge_m = triangle_m.getOppositeNode(edge_m.n1(), edge_m.n2());

            // Triangle current is from triangle with smaller index to larger index
            // First triangle is positive, second triangle is negative (they are sorted)
            double sign_m = RWGBasisFunctionSign(tmIndex, boundingTrianglesIndexList_m.at((triangle_mIndex + 1) % boundingTrianglesIndexList_m.size()));
//            double sign_m = -(double)triangle_mIndex*2.0 + 1.0;
            assert((sign_m == -1.0) || (sign_m == 1.0));

            std::vector<Quadrature::WeightedPoint> weightedPointsTest  = quadrature.getTriangleGaussianQuadraturePoints(triangle_m, m_numberOfTestIntegrationPoints);

            Quadrature::WeightedPoint wpm;
            Quadrature::WeightedPoint wpn;
            std::vector<Quadrature::WeightedPoint> weightedPointsSource;
            arma::cx_vec A (3);
            std::complex<double> B {0.0, 0.0};
            std::complex<double> G0weight;
            Node::vec rho_n_qp;
            Node::vec rho_m_qp;
            Triangle triangle_n;
            Node nodeOppositeToEdge_n;
            Node r_m_qp;
            Node r_n_qp;
            // Loop over qudrature points (to perform outer integration)
            for (unsigned qpmIndex = 0; qpmIndex < weightedPointsTest.size() ; ++qpmIndex)
            {
                wpm = weightedPointsTest.at(qpmIndex);

                // Vector to point from the origin
                r_m_qp = wpm.node;
                rho_m_qp = sign_m * (r_m_qp.getVec() - nodeOppositeToEdge_m.getVec());

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
                        auto tnIndex = boundingTrianglesIndexList_n.at(triangle_nIndex);
                        triangle_n = m_tContainer.at(tnIndex);

                        //set triangle n1 to be the vertex opposite to the edge
                        nodeOppositeToEdge_n = triangle_n.getOppositeNode(edge_n.n1(), edge_n.n2());

                        // First triangle is positive, second triangle is negative
                        double sign_n = RWGBasisFunctionSign(tnIndex, boundingTrianglesIndexList_n.at((triangle_nIndex + 1) % boundingTrianglesIndexList_n.size()));
                        assert((sign_n == -1.0) || (sign_n == 1.0));

                        if (triangle_n.centre().distance(r_m_qp) <  accurateIntegrationDistance)
                        {
                            weightedPointsSource = quadrature.RAR1S(triangle_n, r_m_qp, m_numberOfSourceIntegrationPoints);
                        }
                        else
                        {
                            weightedPointsSource = quadrature.getTriangleGaussianQuadraturePoints(triangle_n, m_numberOfSourceIntegrationPoints);
                        }

                        A.fill({0, 0});
                        B = 0.0;

                        // Loop over qudrature points (to perform integration)
                        for (unsigned qpnIndex = 0; qpnIndex < weightedPointsSource.size() ; ++qpnIndex)
                        {
                            wpn = weightedPointsSource.at(qpnIndex);

                            // Vector to point from the origin
                            r_n_qp = wpn.node;
                            rho_n_qp = sign_n * (r_n_qp.getVec() - nodeOppositeToEdge_n.getVec());

                            double R = r_m_qp.distance(r_n_qp);
                            G0weight = G0(R) * wpn.weight;

                            // CRC: Should use real names for quantities (vector potential and scalar potential)
                            // Calculate A
                            A += rho_n_qp * G0weight;

                            // Calculate B
                            B += G0weight;
                        }
                        std::complex<double> Adotfm = RWGBasisFunction(triangle_m, edge_m) * RWGBasisFunction(triangle_n, edge_n) * dot(A, rho_m_qp) * wpm.weight;
                        std::complex<double> tmp = divRWGBasisFunction(triangle_m, edge_m, sign_m) * divRWGBasisFunction(triangle_n, edge_n, sign_n) * B * wpm.weight;

                        Zmatrix(m, n) += m_k*Adotfm - tmp/m_k;
                    }
                }
            }
        }
    }
    Zmatrix *= EMconst::j * EMconst::Z0;

    return Zmatrix;
}

arma::cx_vec MoM::calculateRHS(PlaneWave pw)
{
    Quadrature quadrature;

    EdgeContainer eContainer(m_tContainer);
    eContainer.buildNonboundaryEdgeList();

    arma::cx_vec Vvector((unsigned)eContainer.size());
    Vvector.fill({0.0, 0.0});

    // For each edge/basis function m
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        auto edge_m = eContainer.at(m);
        auto triangleList_m = edge_m.getTriangles();

        Node::vec rho_m_qp;
        Node r_m_qp;
        // Each triangle at the edge (edge always only has two triangles - we don't support junctions)
        for (unsigned tmindex = 0; tmindex < 2 ; ++tmindex)
        {
            Triangle triangle_m = m_tContainer.at(triangleList_m.at(tmindex));
            triangle_m.setOppositeEdge(eContainer.at(m).n1(), eContainer.at(m).n2()); // Set the remaining vertex as the first node of the triangle

            double sign_m = -(double)tmindex*2.0 + 1.0; // First triangle is positive, second triangle is negative
            assert((sign_m == -1.0) || (sign_m == 1.0));

            std::vector<Quadrature::WeightedPoint> weightedPointsTest = quadrature.getTriangleGaussianQuadraturePoints(triangle_m, m_numberOfSourceIntegrationPoints);

            for (unsigned qpmIndex = 0; qpmIndex < weightedPointsTest.size() ; ++qpmIndex)
            {
                Quadrature::WeightedPoint wpm = weightedPointsTest.at(qpmIndex);

                // Vector to point from the origin
                r_m_qp = wpm.node;
                rho_m_qp = sign_m * (r_m_qp.getVec() - triangle_m.n1().getVec());

                pw.setFieldPoint(r_m_qp);
                Vvector(m) += RWGBasisFunction(triangle_m, edge_m) * arma::dot(pw.getElectricField().getVec(), rho_m_qp) * wpm.weight;
           }
        }
    }

    return Vvector;
}


void MoM::writeCurrentsToOS(std::string fname, arma::cx_vec solutionMatrix) const
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
            std::vector<EdgeContainer::SizeType> edgeIndices = eContainer.getEdgeIndecesOnTriangle(tIndex);

            Triangle T = m_tContainer.at(tIndex);
//            std::cout << std::endl << "Triangle " << tIndex+1 << ": " << std::endl;
//            T.n1().print();
//            T.n2().print();
//            T.n3().print();

            for (auto nIndex=0; nIndex < 3 ; ++nIndex)
            {
                Node n;
                n = T[nIndex];
//                std::cout << "  n" << nIndex << ": " << n.x() << " , " << n.y() << " , " << n.z() << std::endl;

                for (EdgeContainer::SizeType eIndex=0; eIndex < edgeIndices.size() ; ++eIndex)
                {
                    // If edge contains the node of interest
                    if (eContainer.at(edgeIndices.at(eIndex)).n1() == n ||
                        eContainer.at(edgeIndices.at(eIndex)).n2() == n )
                    {
                        Edge e = eContainer.at(edgeIndices.at(eIndex));
                        std::vector<TriangleContainer::SizeType> triangleIndeces = e.getTriangles();

                        TriangleContainer::SizeType otherTriangleIndex = triangleIndeces.at(0);
                        if (tIndex == otherTriangleIndex)
                        {
                            otherTriangleIndex = triangleIndeces.at(1);
                        }
                        double sign = RWGBasisFunctionSign(tIndex, otherTriangleIndex);

                        Triangle Tl = T;
                        Tl.setOppositeEdge(e.n1(), e.n2());
                        Node nodeOppositeToEdge = Tl.n1();

                        Node rho = n - nodeOppositeToEdge;
                        Node current = sign*rho*RWGBasisFunction(Tl, e);

                        NearFieldValue val;
                        val = nodeCurrents[nIndex];
                        val.setX(val.getX() + current.x()*solutionMatrix(edgeIndices.at(eIndex)));
                        val.setY(val.getY() + current.y()*solutionMatrix(edgeIndices.at(eIndex)));
                        val.setZ(val.getZ() + current.z()*solutionMatrix(edgeIndices.at(eIndex)));
//                        std::cout << val.getX() << "  " << val.getY() << "   " << val.getZ() << std::endl;
                        nodeCurrents[nIndex] = val;
                    }
                }
            }

            // Currents at the centre
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

            // Write out current at the node
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
            file << std::endl;
        }
        file << std::endl;
        file.close();
    }
}

