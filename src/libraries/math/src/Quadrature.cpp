#include "Quadrature.hpp"
#include "ComplexMatrix.hpp"
#include "Edge.hpp"
#include "assert.h"
#include <iostream>

Quadrature::Quadrature()
{
    //ctor
}

Quadrature::~Quadrature()
{
    //dtor
}

Quadrature::WeightedPoint1DList_type Quadrature::get1DGaussianQuadraturePoints(unsigned maxNumberOfPoints)
{

    assert(maxNumberOfPoints > 0);

    std::vector<std::vector<double>> pointsWithWeights;
    std::vector<unsigned> validNumber {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    unsigned numberOfIntegrationPoints;

    std::vector<unsigned>::const_iterator validPointIt;
    for (validPointIt = validNumber.cbegin(); validPointIt < validNumber.cend(); ++validPointIt)
    {
        if ( *validPointIt <= maxNumberOfPoints )
        {
            numberOfIntegrationPoints = *validPointIt;
        }
    }

    switch (numberOfIntegrationPoints)
    {
        case 1:
            pointsWithWeights.push_back({0.0, 2.0});
            break;
        case 2:
            pointsWithWeights.push_back({-0.5773502691896257, 1.0});
            pointsWithWeights.push_back({ 0.5773502691896257, 1.0});
            break;
        case 3:
            pointsWithWeights.push_back({ 0.0, 0.8888888888888888});
            pointsWithWeights.push_back({-0.7745966692414834, 0.5555555555555556});
            pointsWithWeights.push_back({ 0.7745966692414834, 0.5555555555555556});
            break;
        case 4:
            pointsWithWeights.push_back({-0.3399810435848563, 0.6521451548625461});
            pointsWithWeights.push_back({ 0.3399810435848563, 0.6521451548625461});
            pointsWithWeights.push_back({-0.8611363115940526, 0.3478548451374538});
            pointsWithWeights.push_back({ 0.8611363115940526, 0.3478548451374538});
            break;
        case 5:
            pointsWithWeights.push_back({ 0.0, 0.5688888888888889});
            pointsWithWeights.push_back({-0.5384693101056831, 0.4786286704993665});
            pointsWithWeights.push_back({ 0.5384693101056831, 0.4786286704993665});
            pointsWithWeights.push_back({-0.9061798459386640, 0.2369268850561891});
            pointsWithWeights.push_back({ 0.9061798459386640, 0.2369268850561891});
            break;
        case 6:
            pointsWithWeights.push_back({ 0.6612093864662645, 0.3607615730481386});
            pointsWithWeights.push_back({-0.6612093864662645, 0.3607615730481386});
            pointsWithWeights.push_back({-0.2386191860831969, 0.4679139345726910});
            pointsWithWeights.push_back({ 0.2386191860831969, 0.4679139345726910});
            pointsWithWeights.push_back({-0.9324695142031521, 0.1713244923791704});
            pointsWithWeights.push_back({ 0.9324695142031521, 0.1713244923791704});
            break;
        case 7:
            pointsWithWeights.push_back({0.0, 0.4179591836734694});
            pointsWithWeights.push_back({ 0.4058451513773972, 0.3818300505051189});
            pointsWithWeights.push_back({-0.4058451513773972, 0.3818300505051189});
            pointsWithWeights.push_back({ 0.7415311855993945, 0.2797053914892766});
            pointsWithWeights.push_back({-0.7415311855993945, 0.2797053914892766});
            pointsWithWeights.push_back({ 0.9491079123427585, 0.1294849661688697});
            pointsWithWeights.push_back({-0.9491079123427585, 0.1294849661688697});
            break;
        case 8:
            pointsWithWeights.push_back({-0.96028985649754, 0.10122853629038});
            pointsWithWeights.push_back({-0.79666647741362, 0.22238103445338});
            pointsWithWeights.push_back({-0.52553240991633, 0.31370664587788});
            pointsWithWeights.push_back({-0.18343464249565, 0.36268378337836});
            pointsWithWeights.push_back({ 0.18343464249565, 0.36268378337836});
            pointsWithWeights.push_back({ 0.52553240991633, 0.31370664587788});
            pointsWithWeights.push_back({ 0.79666647741362, 0.22238103445338});
            pointsWithWeights.push_back({ 0.96028985649754, 0.10122853629038});
            break;
        case 9:
            pointsWithWeights.push_back({-0.96816023950763, 0.08127438836157});
            pointsWithWeights.push_back({-0.83603110732663, 0.18064816069487});
            pointsWithWeights.push_back({-0.61337143270059, 0.26061069640291});
            pointsWithWeights.push_back({-0.32425342340381, 0.31234707704002});
            pointsWithWeights.push_back({ 0.0, 0.33023935500125});
            pointsWithWeights.push_back({ 0.32425342340381, 0.31234707704001});
            pointsWithWeights.push_back({ 0.61337143270059, 0.26061069640292});
            pointsWithWeights.push_back({ 0.83603110732663, 0.18064816069487});
            pointsWithWeights.push_back({ 0.96816023950763, 0.08127438836157});
            break;
        case 10:
        default:
            pointsWithWeights.push_back({-0.97390652851717, 0.06667134430869});
            pointsWithWeights.push_back({-0.86506336668899, 0.14945134915058});
            pointsWithWeights.push_back({-0.67940956829902, 0.21908636251599});
            pointsWithWeights.push_back({-0.43339539412925, 0.26926671930998});
            pointsWithWeights.push_back({-0.14887433898163, 0.29552422471476});
            pointsWithWeights.push_back({ 0.14887433898163, 0.29552422471475});
            pointsWithWeights.push_back({ 0.43339539412925, 0.26926671930999});
            pointsWithWeights.push_back({ 0.67940956829902, 0.21908636251599});
            pointsWithWeights.push_back({ 0.86506336668899, 0.14945134915058});
            pointsWithWeights.push_back({ 0.97390652851717, 0.06667134430869});
            break;
    }

    WeightedPoint1DList_type weightedPoints;
    WeightedPoint1D wp;

    std::vector<std::vector<double>>::const_iterator it;
    for (it = pointsWithWeights.begin(); it < pointsWithWeights.end(); ++it)
    {
        wp.node =   it->at(0);
        wp.weight = it->at(1);
        weightedPoints.push_back(wp);
    }

    return weightedPoints;
}

Quadrature::WeightedPointList_type Quadrature::getTriangleSimplexGaussianQuadraturePoints(unsigned maxNumberOfPoints)
{

    assert(maxNumberOfPoints > 0);

    std::vector<std::vector<double>> pointsWithWeights;
    std::vector<unsigned> validNumber {1, 3, 4, 6, 7, 12, 13, 16, 25, 33};
    unsigned numberOfIntegrationPoints;

    std::vector<unsigned>::const_iterator validPointIt;
    for (validPointIt = validNumber.cbegin(); validPointIt < validNumber.cend(); ++validPointIt)
    {
        if ( *validPointIt <= maxNumberOfPoints )
        {
            numberOfIntegrationPoints = *validPointIt;
        }
    }

    switch (numberOfIntegrationPoints)
    {
        case 1:
            pointsWithWeights.push_back({1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0});
            break;
        case 3:
            pointsWithWeights.push_back({2.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0});
            break;
        case 4:
            pointsWithWeights.push_back({1.0/3, 1.0/3, 1.0/3, -0.5625});
            pointsWithWeights.push_back({0.6, 0.2, 0.2, 0.520833333333333});
            break;
        case 6:
            pointsWithWeights.push_back({0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322});
            pointsWithWeights.push_back({0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011});
            break;
        case 7:
            pointsWithWeights.push_back({1.0/3, 1.0/3, 1.0/3, 0.225});
            pointsWithWeights.push_back({0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827});
            pointsWithWeights.push_back({0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506});
            break;
        case 12:
            pointsWithWeights.push_back({0.873821971016996, 0.063089014491502, 0.063089014491502, 0.050844906370207});
            pointsWithWeights.push_back({0.501426509658179, 0.249286745170910, 0.249286745170910, 0.116786275726379});
            pointsWithWeights.push_back({0.636502499121399, 0.310352451033785, 0.053145049844816, 0.082851075618374});
            pointsWithWeights.push_back({0.636502499121399, 0.053145049844816, 0.310352451033785, 0.082851075618374});
            break;
        case 13:
            pointsWithWeights.push_back({1.0/3, 1.0/3, 1.0/3, -0.149570044467682});
            pointsWithWeights.push_back({0.479308067841920, 0.260345966079040, 0.260345966079040, 0.175615257433208});
            pointsWithWeights.push_back({0.869739794195568, 0.065130102902216, 0.065130102902216, 0.053347235608838});
            pointsWithWeights.push_back({0.048690315425316, 0.312865496004874, 0.638444188569810, 0.077113760890257});
            pointsWithWeights.push_back({0.048690315425316, 0.638444188569810, 0.312865496004874, 0.077113760890257});
            break;
        case 16:
            pointsWithWeights.push_back({1.0/3, 1.0/3, 1.0/3, 0.144315607677787});
            pointsWithWeights.push_back({0.081414823414554, 0.459292588292723, 0.459292588292723, 0.095091634267284});
            pointsWithWeights.push_back({0.658861384496478, 0.170569307751761, 0.170569307751761, 0.103217370534718});
            pointsWithWeights.push_back({0.898905543365938, 0.050547228317031, 0.050547228317031, 0.032458497623198});
            pointsWithWeights.push_back({0.008394777409958, 0.263112829634638, 0.728492392955404, 0.027230314174435});
            pointsWithWeights.push_back({0.008394777409958, 0.728492392955404, 0.263112829634638, 0.027230314174435});
            break;
        case 25:
            pointsWithWeights.push_back({1.0/3, 1.0/3, 1.0/3, 0.090817990382754});
            pointsWithWeights.push_back({0.028844733232685, 0.485577633383657, 0.485577633383657, 0.036725957756467});
            pointsWithWeights.push_back({0.781036849029926, 0.109481575485037, 0.109481575485037, 0.045321059435528});
            pointsWithWeights.push_back({0.141707219414880, 0.307939838764121, 0.550352941820999, 0.072757916845420});
            pointsWithWeights.push_back({0.141707219414880, 0.550352941820999, 0.307939838764121, 0.072757916845420});
            pointsWithWeights.push_back({0.025003534762686, 0.246672560639903, 0.728323904597411, 0.028327242531057});
            pointsWithWeights.push_back({0.025003534762686, 0.728323904597411, 0.246672560639903, 0.028327242531057});
            pointsWithWeights.push_back({0.009540815400299, 0.066803251012200, 0.923655933587500, 0.009421666963733});
            pointsWithWeights.push_back({0.009540815400299, 0.923655933587500, 0.066803251012200, 0.009421666963733});
            break;
        case 33:
        default:
            pointsWithWeights.push_back({0.023565220452390, 0.488217389773805, 0.488217389773805, 0.025731066440455});
            pointsWithWeights.push_back({0.120551215411079, 0.439724392294460, 0.439724392294460, 0.043692544538038});
            pointsWithWeights.push_back({0.457579229975768, 0.271210385012116, 0.271210385012116, 0.062858224217885});
            pointsWithWeights.push_back({0.744847708916828, 0.127576145541586, 0.127576145541586, 0.034796112930709});
            pointsWithWeights.push_back({0.957365299093579, 0.021317350453210, 0.021317350453210, 0.006166261051559});
            pointsWithWeights.push_back({0.115343494534698, 0.608943235779788, 0.275713269685514, 0.040371557766381});
            pointsWithWeights.push_back({0.115343494534698, 0.275713269685514, 0.608943235779788, 0.040371557766381});
            pointsWithWeights.push_back({0.022838332222257, 0.695836086787803, 0.281325580989940, 0.022356773202303});
            pointsWithWeights.push_back({0.022838332222257, 0.281325580989940, 0.695836086787803, 0.022356773202303});
            pointsWithWeights.push_back({0.025734050548330, 0.858014033544073, 0.116251915907597, 0.017316231108659});
            pointsWithWeights.push_back({0.025734050548330, 0.116251915907597, 0.858014033544073, 0.017316231108659});
            break;
    }

    WeightedPointList_type weightedPoints;
    WeightedPoint wp;

    int iter = 1;

    std::vector<std::vector<double>>::const_iterator it;
    for (it = pointsWithWeights.begin(); it < pointsWithWeights.end(); ++it)
    {
        if (it->at(0) != it->at(1))
        {
            iter = 3;
        }
        else
        {
            iter = 1;
        }

        for (int j=0; j<iter; ++j)
        {
            wp.node.set(it->at((0+j)%3),
                        it->at((1+j)%3),
                        it->at((2+j)%3));
            wp.weight = it->at(3);
            weightedPoints.push_back(wp);
        }
    }

    return weightedPoints;
}

/* Functions for Radial Angular R1 Sqrt */
double Quadrature::RAR1SwFromXY (double x, double y)
{
    return std::asinh(x/y);
//    return std::log(std::tan(atan2(y, x)/2.0));
}
double Quadrature::RAR1SqFromYWZ (double y, double w, double z)
{
    return std::sqrt(std::sqrt(std::pow(std::cosh(w),2)*y*y + z*z) - std::fabs(z));
}
double Quadrature::RAR1SyFromQWZ (double q, double w, double z)
{
    return std::sqrt((std::pow((q*q + std::fabs(z)),2) - z*z)/std::pow(std::cosh(w),2));
}
double Quadrature::RAR1SxFromYW (double y, double w)
{
    return y*std::sinh(w);
}

ComplexMatrix Quadrature::getXrotationMatrix(double theta)
{
    ComplexMatrix rotation {3}; // This should be a double matrix and not a complex matrix.
    rotation(0,0) =  1.0;
    rotation(1,1) =  cos(theta);
    rotation(1,2) = -sin(theta);
    rotation(2,1) =  sin(theta);
    rotation(2,2) =  cos(theta);
    return rotation;
}

ComplexMatrix Quadrature::getZrotationMatrix(double phi)
{
    ComplexMatrix rotation {3}; // This should be a double matrix and not a complex matrix.
    rotation(0,0) =  cos(phi);
    rotation(0,1) = -sin(phi);
    rotation(1,0) =  sin(phi);
    rotation(1,1) =  cos(phi);
    rotation(2,2) =  1.0;
    return rotation;
}

Quadrature::WeightedPointList_type Quadrature::RAR1S_2D(Triangle T, double cruxZoffset, unsigned maxNumberOfPoints)
{

    assert(maxNumberOfPoints > 0);
    assert(T.n1().z() == 0);
    assert(T.n2().z() == 0);
    assert(T.n1().z() == 0);

    WeightedPointList_type weightedPoints;
    WeightedPoint wp;
    WeightedPoint1DList_type weightedPoints1D;

    weightedPoints1D = get1DGaussianQuadraturePoints(maxNumberOfPoints);

    Node tempV1 = T.n1() - T.n3();
    Node tempV2 = T.n2() - T.n3();
    Node pVec   = tempV2 - tempV1;

    double phi = atan2(-pVec.y(), pVec.x());
    ComplexMatrix Zrotation = getZrotationMatrix(phi);
    ComplexMatrix invZrotation = Zrotation.inverse();

    tempV1 = tempV1.transform(Zrotation);
    tempV2 = tempV2.transform(Zrotation);

    // Flip triangle if upside down
    double flip = 1.0;
    if (tempV1.y() < 0.0)
    {
        flip = -1.0;
    }
    tempV1.set(tempV1.x(), tempV1.y()*flip, tempV1.z());
    tempV2.set(tempV2.x(), tempV2.y()*flip, tempV2.z());

    if (tempV1.x() > tempV2.x())
    {
        Node tempNode = tempV1;
        tempV1 = tempV2;
        tempV2 = tempNode;
    }

    // The triangle now has node 3 at the origin and nodes 1 and 2 have positive y-values.
    // Node 1 is the left most node and node 2 is the right most node.

    double wLower = RAR1SwFromXY(tempV1.x(), tempV1.y());
    double wUpper = RAR1SwFromXY(tempV2.x(), tempV2.y());
    double wRange = wUpper - wLower;

    for (auto itii = weightedPoints1D.begin(); itii < weightedPoints1D.end(); ++itii)
    {
        double wPoint = wLower + 0.5*wRange*(1.0 + itii->node);

        for (auto itjj = weightedPoints1D.begin(); itjj < weightedPoints1D.end(); ++itjj)
        {
            double qLower = RAR1SqFromYWZ(0.0, wPoint, std::fabs(cruxZoffset));
            double qUpper = RAR1SqFromYWZ(tempV1.y(), wPoint, cruxZoffset);
            double qRange = qUpper - qLower;

            double qPoint = qLower + 0.5*qRange*(1.0 + itjj->node);
            double yPoint = RAR1SyFromQWZ(qPoint, wPoint, cruxZoffset);
            double xPoint = RAR1SxFromYW(yPoint, wPoint);

            Node tempNewPoint {xPoint, flip*yPoint, 0.0};
            tempNewPoint = tempNewPoint.transform(invZrotation);
            tempNewPoint += T.n3();
            wp.node   = tempNewPoint;
            wp.weight = itii->weight * itjj->weight * qRange * wRange * 0.25 * 2.0 * qPoint / cosh(wPoint);
            std::cout << wp.weight << "  " << itii->weight << "  " << itjj->weight << "  " << qPoint << "  " << wPoint << std::endl;
            std::cout << qRange << "  " << wRange << "  " << qPoint << std::endl;
            weightedPoints.push_back(wp);
        }
    }

    return weightedPoints;
}

Quadrature::WeightedPointList_type Quadrature::RAR1S(const Triangle& T, Node observationPoint, unsigned maxNumberOfPoints)
{
    WeightedPointList_type weightedPoints;
    WeightedPointList_type weightedPointsSingleTriangle;
    WeightedPoint wp;

    Node triNormal = T.normal();
    if (triNormal.z() < 0.0)
    {
        triNormal *= 1.0;
    }

    double theta = acos(triNormal.z());
    double phi   = atan2(triNormal.y(), triNormal.x());

    ComplexMatrix Zrotation = getZrotationMatrix(phi);
    ComplexMatrix Xrotation = getXrotationMatrix(theta);
    ComplexMatrix rotationMatrix = Xrotation*Zrotation;
    ComplexMatrix invRotationMatrix = rotationMatrix.inverse();

    Triangle newT = T.transform(rotationMatrix);
    double zOffset = newT.n1().z();
    newT.n1().set(newT.n1().x(), newT.n1().y(), 0.0);
    newT.n2().set(newT.n2().x(), newT.n2().y(), 0.0);
    newT.n3().set(newT.n3().x(), newT.n3().y(), 0.0);

    Node newObservationPoint = observationPoint.transform(rotationMatrix);
    Node observationZOffsetNode {0.0, 0.0, newObservationPoint.z()};
    Node newObservationPointNoZ = newObservationPoint - observationZOffsetNode;

    Node newTriCentre = newT.centre();

    // Split the triangle at the projection of the observation point
    for (unsigned tIndex = 0; tIndex < 3 ; ++tIndex)
    {
        Node edge     = newT[(tIndex+1)%3] - newT[tIndex];
        Node toCentre = newTriCentre - newT[tIndex];
        Node toObs    = newObservationPointNoZ - newT[tIndex];

        double sign = edge.cross(toCentre).z() * edge.cross(toObs).z();

        // RAR1S_2D
        weightedPointsSingleTriangle = RAR1S_2D(newT, zOffset, maxNumberOfPoints);

        // For each point retruned by RAR1S_2D, add actual point
        for (auto it = weightedPointsSingleTriangle.begin(); it < weightedPointsSingleTriangle.end(); ++it)
        {
            wp = *it;
            wp.node.set(wp.node.x(), wp.node.y(), wp.node.z() + zOffset);
            double R = newObservationPoint.distance(wp.node);
            wp.node = wp.node.transform(invRotationMatrix);
            wp.weight = wp.weight * R * sign; // Not sure how / why this is done or if it is correct

            weightedPoints.push_back(wp);
        }
    }

    return weightedPoints;
}
