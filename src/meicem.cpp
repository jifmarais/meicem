#include "Node.hpp"

#include "NodeContainer.hpp"
#include "ComplexMatrix.hpp"
#include "TriangleContainer.hpp"
#include "reader_wilco_input.hpp"
#include "MoM.hpp"
//#include "math.h"
#include "assert.h"
#include <iostream>
//#include <iterator>
//#include <complex>
//#include <cmath>

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
    reader.setFile(baseTestFilesDirectory + "input_cononical_2basisfunction.txt");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    double freq = reader.getFrequency();

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    ComplexMatrix Zmatrix = MoMSetup.fillZmatrixTriangle();
    ComplexMatrix Vvector = MoMSetup.calculateRHS();

    std::cout << std::endl << "Zmatrix" << std::endl;
    Zmatrix.print();

    std::cout << std::endl << "Vvector" << std::endl;
    Vvector.print();

//    (Zmatrix.inverse()*Zmatrix).print();
    std::cout << std::endl << "Coefficients" << std::endl;
    ComplexMatrix Vsolution {(unsigned)Vvector.getRowCount(), 1};
    Vsolution = Zmatrix.inverse()*Vvector;
    Vsolution.print();

    // Calculate currents (and export them)

    // Calculate near fields and export them

    return EXIT_SUCCESS;
}

