#include "Node.hpp"
#include "NodeContainer.hpp"
#include "Vector.hpp"
#include "TriangleContainer.hpp"
#include "reader_wilco_input.hpp"
#include <iostream>
#include <iterator>

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

    for (TriangleContainer::SizeType observationIndex = 0; observationIndex < tContainer.size() ; ++observationIndex)
    {
        Triangle tObserver {tContainer.at(observationIndex)};

        for (TriangleContainer::SizeType sourceIndex = 0; sourceIndex < tContainer.size() ; ++sourceIndex)
        {
            Triangle tSource {tContainer.at(sourceIndex)};
            std::cout << "Distance ("
                      << observationIndex
                      << ", " << sourceIndex
                      << "):" << std::to_string( Node::distance(tObserver.centre(), tSource.centre()) * 1e3 )
                      << "    " << tContainer.hasCommonNode(observationIndex, sourceIndex)
                      << std::endl;
        }

    }

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
