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
	, propNumberOfStreamlines("nrOfStreamlines", "Number of Streamlines")

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

	addProperty(propNumberOfStreamlines);

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

dvec2 StreamlineIntegrator::drawStreamline(std::shared_ptr<BasicMesh> &mesh,
                                          VectorField2 vectorField,
                                          std::vector<BasicMesh::Vertex> &vertices, dvec2 seedPoint,
                                          float stepsize) {
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);

    indexBufferRK->add(static_cast<std::uint32_t>(0));

    for (int i = 0; i < propMaxIntegrationSteps; i++) {
        if (propDirectionField == 1) {
        }
        if (propStepsize * i > propMaxArcLength) {
            break;
        }
        seedPoint = Integrator::RK4(vectorField, seedPoint, stepsize);
        indexBufferRK->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back(
            {vec3(seedPoint[0], seedPoint[1], 0), vec3(1), vec3(1), vec4(1, 0, 0, 1)});
	}

	return seedPoint;
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

	auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
	indexBufferPoints->add(static_cast<std::uint32_t>(0));
    if (propSeedMode.get() == 0) {
        // Draw start point
        vec2 s_pos = propStartPoint.get();
        vertices.push_back({vec3(s_pos.x, s_pos.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});

        // TODO: Create one stream line from the given start point

		//c
		if (propDirectionField == 1) {

			auto normalizedVectorField = VectorField2::createFieldFromVolume(vol);
			
			vec2 currentVector = vectorField.interpolate(s_pos);

			float length = sqrt(currentVector[0] * currentVector[0] + currentVector[1] * currentVector[1]);
			LogProcessorInfo("length: " << length)
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
                vec2 vector = drawStreamline(mesh, vectorField, vertices, s_pos, propIntegration ? -propStepsize : propStepsize);

				//g
				if (vector.x == 0 && vector.y == 0)
				{
					break;
				}

				stepCounter++;
			}
		}
    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
		// !!THIS IS NOT DONE YET!!
		std::vector<BasicMesh::Vertex> random_seeds;
		float width = vectorField.getNumVerticesPerDim()[0] * vectorField.getCellSize()[0] ;
		float height = vectorField.getNumVerticesPerDim()[1] * vectorField.getCellSize()[1];
		for (int i = 0; i < propNumberOfStreamlines; i++)
		{
			float randXDir = ((rand() % 100) / (100.0f)) > 0.5 ? 1 : -1;
			float randX = ((rand() % 100) / (100.0f)) * width * randXDir;
			float randYDir = ((rand() % 100) / (100.0f)) > 0.5 ? 1 : -1;
			float randY = ((rand() % 100) / (100.0f)) * height * randYDir;

			LogProcessorInfo("RAND VEC: " << vec2(randX, randY));
			
			vertices.push_back({ vec3(randX, randY, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1) });

			drawStreamline(mesh, vectorField, vertices, vec2(randX, randY), propStepsize);
		}
		
		// (TODO: Bonus, sample randomly according to magnitude of the vector field)
    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

}  // namespace inviwo
