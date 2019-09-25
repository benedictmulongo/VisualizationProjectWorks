/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/util/glm.h>
#include <math.h>

namespace inviwo {

// TODO: Implement a single integration step here
dvec2 Integrator::normalize(const VectorField2& vectorField, const dvec2& pos)
{
    dvec2 new_vec = vectorField.interpolate(pos);
    double length = sqrt(new_vec[0]*new_vec[0] + new_vec[1]* new_vec[1] );
    length = (length == 0 ) ? 1 : (1 / length );

    return length * new_vec;
}

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position,
                        const double& stepSize) {
    // Access the vector field with vectorField.interpolate(...)
    dvec2 oneStep = position + stepSize * vectorField.interpolate(position);
    return oneStep;
}

dvec2 Integrator::Euler_stream(const VectorField2& vectorField, const dvec2& position,
                        const double& stepSize) {
    // Access the vector field with vectorField.interpolate(...)
    dvec2 oneStep = position + stepSize * Integrator::normalize(vectorField,position);
    return oneStep;
}

dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position,
                      const double& stepSize) {
    dvec2 k1 = stepSize * vectorField.interpolate(position);
    dvec2 k2 = stepSize * vectorField.interpolate(position + 0.5 * k1);
    dvec2 k3 = stepSize * vectorField.interpolate(position + 0.5 * k2);
    dvec2 k4 = stepSize * vectorField.interpolate(position + k3);

    dvec2 B = 2 * k2;
    dvec2 C = 2 * k3;
    dvec2 sum = k1 + B + C + k4;
    dvec2 oneStep = position + 0.16666666666 * sum;

    Integrator::normalize(vectorField, dvec2(1,1));
    return oneStep;
}


dvec2 Integrator::RK4_stream(const VectorField2& vectorField, const dvec2& position,
const double& stepSize) {

    dvec2 k1 = stepSize * Integrator::normalize(vectorField,position);
    dvec2 k2 = stepSize * Integrator::normalize(vectorField,position + 0.5 * k1);
    dvec2 k3 = stepSize * Integrator::normalize(vectorField,position + 0.5 * k2);
    dvec2 k4 = stepSize * Integrator::normalize(vectorField,position + k3);

    dvec2 B = 2 * k2;
    dvec2 C = 2 * k3;
    dvec2 sum = k1 + B + C + k4;
    dvec2 oneStep = position + 0.16666666666 * sum;

    return oneStep;
}

}  // namespace inviwo
