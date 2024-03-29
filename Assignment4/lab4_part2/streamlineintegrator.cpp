/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , outMesh("meshOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event *e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
	// TODO: Initialize additional properties
	// propertyName("propertyIdentifier", "Display Name of the Propery",
	// default value (optional), minimum value (optional), maximum value (optional),
	// increment (optional)); propertyIdentifier cannot have spaces
    , propIntegration("integration", "Integration")
    , propDirectionField("directionField", "Integration in Direction Field")
    , propStepsize("stepsize", "Step Size", 0.1, 0, 1)
    , propMaxIntegrationSteps("maxintegrationsteps", "Maximum Integration Steps", 50, 0, 1000)
    , propMaxArcLength("maxarclength", "Maximum Arc Length of Streamlines", 100, 0, 1000)
    , propStopAtBoundary("stopAtBoundary", "Stop Integration at Boundary")
    , propStopAtZeros("stopAtZeros", "Stop Integration at Zeros of the Vector Field")
    , propMinVelocity("minVelocity", "Minimum Velocity")

{
    // Register Ports
    addPort(inData);
    addPort(outMesh);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
	//b
    addProperty(propStepsize);

    addProperty(propIntegration);
    propIntegration.addOption("forward", "Forward", 0);
    propIntegration.addOption("backward", "Backward", 1);
    
	addProperty(propDirectionField);

    addProperty(propMaxIntegrationSteps);

    addProperty(propMaxArcLength);

    addProperty(propStopAtBoundary);

    addProperty(propStopAtZeros);

    addProperty(propMinVelocity);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart);
            // util::hide(...)
        } else {
            util::hide(propStartPoint, mouseMoveStart);
            // util::show(...)
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event *event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent *>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to range [0,1]^2
    mousePos = mousePos * 2 - vec2(1, 1);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);

    // The start point should be inside the volume (set maximum to the upper
    // right corner)
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propSeedMode.get() == 0) {
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        // Draw start point
        vec2 startPoint = propStartPoint.get();
        vertices.push_back(
            {vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
        indexBufferPoints->add(static_cast<std::uint32_t>(0));

        // TODO: Create one stream line from the given start point
        dvec2 s_pos = startPoint;
        auto indexBufferRK = mesh->addIndexBuffer(DrawType:: Lines, ConnectivityType::Strip);

        indexBufferRK->add(static_cast<std::uint32_t>(0));

		//c
		if (propDirectionField == 1) {

			auto normalizedVectorField = VectorField2::createFieldFromVolume(vol);
			
			vec2 currentVector = vectorField.interpolate(s_pos);

			float length = sqrt(currentVector[0] * currentVector[0] + currentVector[1] * currentVector[1]);
			//LogProcessorInfo("length: " << length)
			vec2 normalizedVector = currentVector / length;

			normalizedVectorField.setValueAtVertex(s_pos, normalizedVector);
		
			LogProcessorInfo("original: " << currentVector)
			LogProcessorInfo("normalized: " << normalizedVector)
			
			vectorField = normalizedVectorField;
		}
		else {
			vectorField = VectorField2::createFieldFromVolume(vol);
		}

		int stepCounter = 0;
		//f
		for (int i = 0; i < vectorField.getNumVerticesPerDim()[0]; i++) {
			for (int j = 0; j < vectorField.getNumVerticesPerDim()[1]; j++) {
				
				//d
				if (stepCounter >= propMaxIntegrationSteps)
				{
					break;
				}
				
				//e
				if (propStepsize * i > propMaxArcLength) {
					break;
				}

				//h
				vec2 currentVector = vectorField.interpolate(s_pos);
				float velocity = sqrt(currentVector[0] * currentVector[0] + currentVector[1] * currentVector[1]);
				if (velocity < propMinVelocity)
				{
					break;
				}

				//a
				s_pos = Integrator::RK4(vectorField, s_pos, propIntegration ? -propStepsize : propStepsize);

				//g
				if (s_pos[0] == 0 && s_pos[1] == 0)
				{
					break;
				}

				stepCounter++;

				indexBufferRK->add(static_cast<std::uint32_t>(vertices.size()));
				vertices.push_back({ vec3(s_pos[0], s_pos[1], 0), vec3(1), vec3(1), vec4(1, 0, 0, 1) });
			}
		}
    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

}  // namespace inviwo
