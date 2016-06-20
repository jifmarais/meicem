#include "Quadrature.hpp"
#include "assert.h"

Quadrature::Quadrature()
{
    //ctor
}

Quadrature::~Quadrature()
{
    //dtor
}

Quadrature::WeightedPointList_type Quadrature::getTriangleSimplexGaussianQuadraturePoints(unsigned maxNumberOfPoints)
{

    assert(maxNumberOfPoints > 0);

    std::vector<std::vector<double>> pointsWithWeights;
    std::vector<WeightedPoint> weightedPoints;
    WeightedPoint wp;
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

