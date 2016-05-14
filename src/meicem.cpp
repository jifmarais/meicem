#include "Point3D.hpp"
#include "Point3DContainer.hpp"
#include <iostream>
#include <iterator>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main( int argc, char *argv[] )
{
    Point3DContainer L1;
    Point3D p1;
    Point3DContainer::SizeType ii;
    Point3DContainer::SizeType count {50};

    // Add a bunch of points
    for ( ii=0; ii < count ; ++ii )
    {
        p1.set(1.0*ii, 2.0*ii, 3.0*ii);

        L1.addPoint(p1);
    }

    // Points should not be added multiple times 
    for ( ii=0; ii < count ; ++ii )
    {
        p1.set(1.0*ii, 2.0*ii, 3.0*ii);

        L1.addPoint(p1);
    }

    std::cout << "Hello world!" << std::endl;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("compression", po::value<double>(), "set compression level")
            ;

        po::variables_map vm;    
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }

        if (vm.count("compression")) {
            std::cout << "Compression level was set to " 
                << vm["compression"].as<double>() << ".\n";
        } else {
            std::cout << "Compression level was not set.\n";
        }
    }   
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }   
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }   

    return EXIT_SUCCESS;
}
