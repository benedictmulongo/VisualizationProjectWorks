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

namespace inviwo {

// TODO: Implement a single integration step here

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position, const double& stepSize)
{
    //Access the vector field with vectorField.interpolate(...)

    dvec2 oneStep = position + stepSize * vectorField.interpolate(position);

    return oneStep;
}

dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const double& stepSize)
{
    dvec2 k1 =  stepSize * vectorField.interpolate(position);
    dvec2 k2 =  stepSize * vectorField.interpolate(position + 0.5*k1);
    dvec2 k3 =  stepSize * vectorField.interpolate(position + 0.5*k2);
    dvec2 k4 =  stepSize * vectorField.interpolate(position + k3);

   // dvec2 oneStep = vec2(5,6);
    
    dvec2 B = 2*k2;
    dvec2 C = 2*k3;
    dvec2 somme = k1 + B + C + k4 ;
    dvec2 oneStep = position + 0.16666666666*somme;

    return oneStep;
}

}  // namespace inviwo
