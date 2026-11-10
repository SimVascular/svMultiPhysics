// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file SymmetricTriangleQuadrature.cpp
 * @brief Fixed, contained Dunavant tables and ordered symmetry expansion.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/SymmetricTriangleQuadrature.h"

#include "FE/Common/FEException.h"

#include <algorithm>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

// One orbit of a symmetric rule: the points whose barycentric coordinates
// (1 - x - y, x, y) are the distinct permutations of one triple, all with the
// same weight. Members follow the columns of Dunavant's Appendix II.
struct Orbit {
    std::size_t size; // 1: centroid; 3: permutations of (a, a, 1 - 2a); 6: of (a, b, 1 - a - b).
    double weight;    // Published weight of each point (a rule's weights sum to 1); halved on expansion.
    double a = 0.0;   // 3 points: beta (= gamma); 6 points: alpha.
    double b = 0.0;   // 6 points: beta.
};

// kDegreeN is rule N of Dunavant (1985), https://doi.org/10.1002/nme.1620210612,
// in the orbit order and 15-decimal precision of its Appendix II. Repeating
// decimals are written as fractions.
constexpr Orbit kDegree1[] = {{1, 1.0}};
constexpr Orbit kDegree2[] = {{3, 1.0 / 3.0, 1.0 / 6.0}};
constexpr Orbit kDegree3[] = {{1, -0.5625}, {3, 25.0 / 48.0, 0.2}};

constexpr Orbit kDegree4[] = {
    {3, 0.223381589678011, 0.445948490915965},
    {3, 0.109951743655322, 0.091576213509771},
};

constexpr Orbit kDegree5[] = {
    {1, 0.225},
    {3, 0.132394152788506, 0.470142064105115},
    {3, 0.125939180544827, 0.101286507323456},
};

constexpr Orbit kDegree6[] = {
    {3, 0.116786275726379, 0.249286745170910},
    {3, 0.050844906370207, 0.063089014491502},
    {6, 0.082851075618374, 0.053145049844817, 0.310352451033784},
};

constexpr Orbit kDegree7[] = {
    {1, -0.149570044467682},
    {3, 0.175615257433208, 0.260345966079040},
    {3, 0.053347235608838, 0.065130102902216},
    {6, 0.077113760890257, 0.048690315425316, 0.312865496004874},
};

constexpr Orbit kDegree8[] = {
    {1, 0.144315607677787},
    {3, 0.095091634267285, 0.459292588292723},
    {3, 0.103217370534718, 0.170569307751760},
    {3, 0.032458497623198, 0.050547228317031},
    {6, 0.027230314174435, 0.008394777409958, 0.263112829634638},
};

constexpr Orbit kDegree9[] = {
    {1, 0.097135796282799},
    {3, 0.031334700227139, 0.489682519198738},
    {3, 0.077827541004774, 0.437089591492937},
    {3, 0.079647738927210, 0.188203535619033},
    {3, 0.025577675658698, 0.044729513394453},
    {6, 0.043283539377289, 0.036838412054736, 0.221962989160766},
};

constexpr Orbit kDegree10[] = {
    {1, 0.090817990382754},
    {3, 0.036725957756467, 0.485577633383657},
    {3, 0.045321059435528, 0.109481575485037},
    {6, 0.072757916845420, 0.141707219414880, 0.307939838764121},
    {6, 0.028327242531057, 0.025003534762686, 0.246672560639903},
    {6, 0.009421666963733, 0.009540815400299, 0.066803251012200},
};

constexpr Orbit kDegree12[] = {
    {3, 0.025731066440455, 0.488217389773805},
    {3, 0.043692544538038, 0.439724392294460},
    {3, 0.062858224217885, 0.271210385012116},
    {3, 0.034796112930709, 0.127576145541586},
    {3, 0.006166261051559, 0.021317350453210},
    {6, 0.040371557766381, 0.115343494534698, 0.275713269685514},
    {6, 0.022356773202303, 0.022838332222257, 0.281325580989940},
    {6, 0.017316231108659, 0.025734050548330, 0.116251915907597},
};

constexpr Orbit kDegree13[] = {
    {1, 0.052520923400802},
    {3, 0.011280145209330, 0.495048184939705},
    {3, 0.031423518362454, 0.468716635109574},
    {3, 0.047072502504194, 0.414521336801277},
    {3, 0.047363586536355, 0.229399572042831},
    {3, 0.031167529045794, 0.114424495196330},
    {3, 0.007975771465074, 0.024811391363459},
    {6, 0.036848402728732, 0.094853828379579, 0.268794997058761},
    {6, 0.017401463303822, 0.018100773278807, 0.291730066734288},
    {6, 0.015521786839045, 0.022233076674090, 0.126357385491669},
};

constexpr Orbit kDegree14[] = {
    {3, 0.021883581369429, 0.488963910362179},
    {3, 0.032788353544125, 0.417644719340454},
    {3, 0.051774104507292, 0.273477528308839},
    {3, 0.042162588736993, 0.177205532412543},
    {3, 0.014433699669777, 0.061799883090873},
    {3, 0.004923403602400, 0.019390961248701},
    {6, 0.024665753212564, 0.057124757403648, 0.172266687821356},
    {6, 0.038571510787061, 0.092916249356972, 0.336861459796345},
    {6, 0.014436308113534, 0.014646950055654, 0.298372882136258},
    {6, 0.005010228838501, 0.001268330932872, 0.118974497696957},
};

constexpr Orbit kDegree17[] = {
    {1, 0.033437199290803},
    {3, 0.005093415440507, 0.497170540556774},
    {3, 0.014670864527638, 0.482176322624625},
    {3, 0.024350878353672, 0.450239969020782},
    {3, 0.031107550868969, 0.400266239377397},
    {3, 0.031257111218620, 0.252141267970953},
    {3, 0.024815654339665, 0.162047004658461},
    {3, 0.014056073070557, 0.075875882260746},
    {3, 0.003194676173779, 0.015654726967822},
    {6, 0.008119655318993, 0.010186928826919, 0.334319867363658},
    {6, 0.026805742283163, 0.135440871671036, 0.292221537796944},
    {6, 0.018459993210822, 0.054423924290583, 0.319574885423190},
    {6, 0.008476868534328, 0.012868560833637, 0.190704224192292},
    {6, 0.018292796770025, 0.067165782413524, 0.180483211648746},
    {6, 0.006665632004165, 0.014663182224828, 0.080711313679564},
};

constexpr Orbit kDegree19[] = {
    {1, 0.032906331388919},
    {3, 0.010330731891272, 0.489609987073006},
    {3, 0.022387247263016, 0.454536892697893},
    {3, 0.030266125869468, 0.401416680649431},
    {3, 0.030490967802198, 0.255551654403098},
    {3, 0.024159212741641, 0.177077942152130},
    {3, 0.016050803586801, 0.110061053227952},
    {3, 0.008084580261784, 0.055528624251840},
    {3, 0.002079362027485, 0.012621863777229},
    {6, 0.003884876904981, 0.003611417848412, 0.395754787356943},
    {6, 0.025574160612022, 0.134466754530780, 0.307929983880436},
    {6, 0.008880903573338, 0.014446025776115, 0.264566948406520},
    {6, 0.016124546761731, 0.046933578838178, 0.358539352205951},
    {6, 0.002491941817491, 0.002861120350567, 0.157807405968595},
    {6, 0.018242840118951, 0.223861424097916, 0.075050596975911},
    {6, 0.010258563736199, 0.034647074816760, 0.142421601113383},
    {6, 0.003799928855302, 0.010161119296278, 0.065494628082938},
};

struct Table {
    int degree;
    std::span<const Orbit> orbits;
};

constexpr Table kTables[] = {
    {1, kDegree1}, {2, kDegree2}, {3, kDegree3}, {4, kDegree4}, {5, kDegree5},
    {6, kDegree6}, {7, kDegree7}, {8, kDegree8}, {9, kDegree9}, {10, kDegree10},
    {12, kDegree12}, {13, kDegree13}, {14, kDegree14}, {17, kDegree17}, {19, kDegree19},
};

// Every valid request has a table at or above it.
static_assert(kTables[std::size(kTables) - 1].degree == max_symmetric_triangle_exactness());

} // namespace

QuadratureRule make_symmetric_triangle_rule(int requested_exactness)
{
    svmp::check<InvalidArgumentException>(
        requested_exactness >= 0 && requested_exactness <= max_symmetric_triangle_exactness(),
        "Symmetric triangle generator: requested_exactness must be in [0, 19]");
    const Table& table = *std::find_if(std::begin(kTables), std::end(kTables),
        [requested_exactness](const Table& candidate) { return candidate.degree >= requested_exactness; });
    std::size_t count = 0;
    for (const Orbit& orbit : table.orbits) {
        count += orbit.size;
    }
    std::vector<QuadPoint> points;
    std::vector<double> weights;
    points.reserve(count);
    weights.reserve(count);
    for (const auto& [size, weight, a, b] : table.orbits) {
        if (size == 1) {
            points.push_back({1.0 / 3.0, 1.0 / 3.0, 0.0});
        } else if (size == 3) {
            const double c = 1.0 - 2.0 * a;
            points.insert(points.end(), {{a, c, 0.0}, {a, a, 0.0}, {c, a, 0.0}});
        } else {
            const double c = 1.0 - a - b;
            points.insert(points.end(),
                {{b, c, 0.0}, {c, b, 0.0}, {a, c, 0.0}, {c, a, 0.0}, {a, b, 0.0}, {b, a, 0.0}});
        }
        weights.insert(weights.end(), size, 0.5 * weight);
    }
    return QuadratureRule(svmp::CellFamily::Triangle, table.degree, std::move(points), std::move(weights));
}

} // namespace svmp::FE::quadrature
