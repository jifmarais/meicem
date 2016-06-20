#include "Node.hpp"

#include "NodeContainer.hpp"
#include "EdgeContainer.hpp"
#include "TriangleContainer.hpp"
#include "reader_wilco_input.hpp"
#include "Quadrature.hpp"
#include "math.h"
#include "assert.h"
#include <iostream>
#include <iterator>
#include <complex>
#include <cmath>

//#include <boost/program_options.hpp>
//namespace po = boost::program_options;


//int main( int argc, char *argv[] )
int main()
{
    std::cout << "Let the MEICEM mayhem begin..." << std::endl;

    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    WilcoInputReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/reader_wilco_input/test/test_files/";
    reader.setFile(baseTestFilesDirectory + "input.txt");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    Node p1 {0, 0, 2};

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

    double freq = reader.getFrequency();
    const std::complex<double> j (0.0, -1.0);
    const double eps0 {8.85418781761e-12};
    const double mu0  {1.25663706143592e-06};
    const double c0 {299792456.2};

    double omega = 2*M_PI*freq;
    double k = (2*M_PI*freq)/c0;

    std::complex<double> Z_mn {0.0, 0.0};

    EdgeContainer eContainer(tContainer);
    eContainer.buildNonboundaryEdgeList();
    std::cout << "Triangle container size: " << tContainer.size() << std::endl;
    std::cout << "Edge container size:     " << eContainer.size() << std::endl;
    for (unsigned m = 0 ; m < eContainer.size() ; ++m)
    {
        auto triangleList_m = eContainer.at(m).getTriangles();

        for (unsigned tmindex = 0; tmindex < 2 ; ++tmindex)
        {
            Triangle triangle_m = tContainer.at(triangleList_m.at(tmindex));
            triangle_m.setOppositeEdge(eContainer.at(m).n1(), eContainer.at(m).n2());
            Node rmc = triangle_m.centre();
            double sign_m = (tmindex <= 0) ? 1.0 : -1.0;
            Node rho_m_c = (rmc - triangle_m.n1()) * sign_m; // Vector from opposite vertex to centre of triangle (for positive triangle)

            for (unsigned n = 0 ; n < eContainer.size() ; ++n)
            {
                auto triangleList_n = eContainer.at(n).getTriangles();
                // For the positive and negative triangle
                for (unsigned tnindex = 0; tnindex < 2 ; ++tnindex)
                {
                    Triangle triangle_n = tContainer.at(triangleList_n.at(tnindex));
                    triangle_n.setOppositeEdge(eContainer.at(n).n1(), eContainer.at(n).n2());

                    // Calculate quadrature points
                    std::vector<Quadrature::WeightedPoint> weightedPoints;
                    weightedPoints = Quadrature::getTriangleSimplexGaussianQuadraturePoints(3);

                    std::complex<double> Phi_mn {0.0, 0.0};
                    std::complex<double> A_mn {0.0, 0.0};

                    // Loop over qudrature points (to perform integration)
                    for (unsigned ii = 0; ii < weightedPoints.size() ; ++ii)
                    {
                        Quadrature::WeightedPoint wp;
                        wp = weightedPoints.at(ii);
                        Node r = triangle_n.fromSimplex(wp.node); // Need to also supply the edge so that it is oriented correctly

                        // Calculate Rm
                        double Rm = Node::distance(rmc,r);

                        // Calculate Phi_mn
                        double sign_n = (tnindex <= 0) ? 1.0 : -1.0;
                        Phi_mn+= -1.0/(4.0*M_PI*j*omega*eps0) *
                                ( sign_n * eContainer.at(n).length()/triangle_n.area()) *
                                (exp(-1.0*j*k*Rm) / Rm) *
                                (wp.weight/triangle_n.area());

                        // Calculate A_mn
                        A_mn += mu0/(4.0*M_PI) *
                                ( eContainer.at(n).length()/( triangle_n.area() * 2.0) ) *
                                wp.node.x() * // This is suppose to be vector
                                (exp(-1.0*j*k*Rm) / Rm) *
                                (wp.weight/triangle_n.area());

                    }

                    // Add contribution to Z_mn
                    Z_mn += eContainer.at(n).length()* (j*omega * (A_mn *(rho_m_c.x())/2.0) + (Phi_mn) );
                }
            }
        }
    }


//    Triangle t1;
//    Node n1 {0, 0, 2};
//    Node n2 {3, 1, 0};
//    Node n3 {2, 1, 5};
//    t1.set(n1, n2, n3);

//    std::vector<Quadrature::WeightedPoint> weightedPoints;
//    weightedPoints = Quadrature::getTriangleSimplexGaussianQuadraturePoints(3);
//    for (unsigned ii = 0; ii < weightedPoints.size() ; ++ii)
//    {
//        Quadrature::WeightedPoint wp;
//        wp = weightedPoints.at(ii);
//        std::cout << ii+1 << ": "
//                  << wp.node.x() << "  "
//                  << wp.node.y() << "  "
//                  << wp.node.z() << "  "
//                  << wp.weight << std::endl;

//        std::cout << t1.fromSimplex(wp.node).x() << "  "
//                  << t1.fromSimplex(wp.node).y() << "  "
//                  << t1.fromSimplex(wp.node).z() << "  "
//                  << std::endl;
//    }

    return EXIT_SUCCESS;
}

// Some code for later - adding command line argument support
//    try {

//        po::options_description desc("Allowed options");
//        desc.add_options()
//            ("help", "produce help message")
//            ("compression", po::value<double>(), "set compression level")
//            ;

//        po::variables_map vm;
//        po::store(po::parse_command_line(argc, argv, desc), vm);
//        po::notify(vm);

//        if (vm.count("help")) {
//            std::cout << desc << "\n";
//            return 0;
//        }

//        if (vm.count("compression")) {
//            std::cout << "Compression level was set to "
//                << vm["compression"].as<double>() << ".\n";
//        } else {
//            std::cout << "Compression level was not set.\n";
//        }
//    }
//    catch(std::exception& e) {
//        std::cerr << "error: " << e.what() << "\n";
//        return 1;
//    }
//    catch(...) {
//        std::cerr << "Exception of unknown type!\n";
//    }
