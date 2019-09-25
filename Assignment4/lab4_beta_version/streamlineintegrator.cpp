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
#include <math.h>
#include <vector> 
#include <random>
#include <cstdlib> // required for srand(), rand().
using namespace std;

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
    , propGridpoints("gridPoint", "Grid Point", vec2(5, 5), vec2(1,1), vec2(50,50), vec2(1))
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
    , propMaxIntegrationSteps("maxintegrationsteps", "Maximum Integration Steps", 20, 0, 500,1)
    , propNumberRamdomPoints("NumberRamdomPoints", "N. Ramdom Points", 20, 1, 500,1)
    , propMaxArcLength("maxarclength", "Maximum Arc Length of Streamlines", 5, 0, 50,0.01)
    , propStopAtBoundary("stopAtBoundary", "Stop Integration at Boundary")
    , propUniformGrid("UniformGrid", "Uniform Grid")
    , propRandomGrid("RandomGrid", "Random Grid")
    , propShowPoints("ShowPoints", "Show Points")
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
    addProperty(propGridpoints);
    

    // TODO: Register additional properties
	//b
    addProperty(propStepsize);

    addProperty(propIntegration);
    propIntegration.addOption("forward", "Forward", 0);
    propIntegration.addOption("backward", "Backward", 1);
    
	addProperty(propDirectionField);

    addProperty(propMaxIntegrationSteps);
    addProperty(propNumberRamdomPoints);

    addProperty(propMaxArcLength);

    addProperty(propStopAtBoundary);
    addProperty(propUniformGrid);
    addProperty(propRandomGrid);
    addProperty(propShowPoints);

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

bool StreamlineIntegrator::insideDomain(const dvec2 bBoxMin, const dvec2 bBoxMax, const dvec2 pos)
{
    double min_x = bBoxMin[0];
    double min_y = bBoxMin[1];
    double max_x =  bBoxMax[0];
    double max_y = bBoxMax[1];

    //return (left <= pos[0]) && (bottom <= pos[1]) && (pos[0] <= right) && (pos[1] <= top);
    return (min_x <= pos[0]) && (min_y <= pos[1]) && (pos[0] <= max_x) && (pos[1] <= max_y);
}

 dvec2 StreamlineIntegrator::generateRandom(const dvec2 bBoxMin, const dvec2 bBoxMax)
 {
    double min_x = bBoxMin[0];
    double min_y = bBoxMin[1];
    double max_x =  bBoxMax[0];
    double max_y = bBoxMax[1];

    default_random_engine generator;
    uniform_real_distribution<double> distribution_x(min_x, max_x);
    double number_x = distribution_x(generator);
    uniform_real_distribution<double> distribution_y(min_y, max_y);
    double number_y = distribution_y(generator);

    return dvec2(number_x, number_y);
 }

double StreamlineIntegrator::velocity(const VectorField2& vectorField, const dvec2& pos)
{
    dvec2 new_vec = vectorField.interpolate(pos);
    double length = sqrt(new_vec[0]*new_vec[0] + new_vec[1]* new_vec[1] );
    return length;
}

vector<double> StreamlineIntegrator::linspace(const double& min, const double& max, const int N)
{
    double h = (max - min) / (N - 1) ;
    vector<double> x_space(N);

    vector<double>::iterator it;
    double val = min;
    for(it = x_space.begin(); it != x_space.end(); ++it )
    {
        *it = val;
        val = val + h;
    }

    return x_space;
}

void StreamlineIntegrator::drawStreamline(std::shared_ptr<BasicMesh> &mesh,
                                          VectorField2 vectorField,
                                          std::vector<BasicMesh::Vertex> &vertices, dvec2 start_pos) 
{
        
        // - bounding box {xmin, ymin} - {xmax, ymax}
        const dvec2 bBoxMin = vectorField.getBBoxMin();
        const dvec2 bBoxMax = vectorField.getBBoxMax();

        //const dvec2 bBoxMin = vectorField.getMinValue();
        //const dvec2 bBoxMax = vectorField.getMaxValue();

        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);

        if(propShowPoints)
        {
            indexBufferPoints->add(static_cast<std::uint32_t>(0));
            vertices.push_back({vec3(start_pos.x, start_pos.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
        }

        indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(start_pos[0], start_pos[1], 0), vec3(1), vec3(1), vec4(0, 0, 0, 1)});
        
        // TODO: Create one stream line from the given start point

        // This is the maximun of arc length Task. 4.2 (e)
        float streamline_length = 0;
        for(int i = 0; streamline_length < propMaxArcLength; ++i, streamline_length += propStepsize)
        {
            // Stop the integration after N stepsize Task 4.2 (d)
            if(i > propMaxIntegrationSteps) break; 

            // Allow integration in the directionField Task 4.2 (c) RK4_stream
            if (propDirectionField == 1) 
            {
                // Allow integration in forward and backward direction Task 4.2 (a,b)
                double stepdirection = (propIntegration == 0) ? propStepsize : -propStepsize;
                start_pos = Integrator::RK4_stream(vectorField, start_pos, stepdirection);

            }else
            {
                // Allow integration in forward and backward direction Task 4.2 (a,b)
                double stepdirection = (propIntegration == 0) ? propStepsize : -propStepsize;
                start_pos = Integrator::RK4(vectorField, start_pos, stepdirection);
            }
            // Stop the integration at the boundary of the domain Task 4.2 (f)
            if(propStopAtBoundary) 
            {
                if(!insideDomain(bBoxMin, bBoxMax, start_pos)) break;
            }
            // Stop the integration when the velocity become slow Task 4.2 (h)
            if(velocity(vectorField, start_pos) < 1.0e-5) break;
            // ------------------------------------------------- HERE
            // Stop the integration at the zeros of the vector field Task 4.2 (g)
            if(propStopAtZeros)
            {
                //if(start_pos[0] < 1.0e-10 && start_pos[1] < 1.0e-10) break;
                // of 
                dvec2 new_vec = vectorField.interpolate(start_pos);
                if(new_vec[0]== 0 && new_vec[1] == 0) break;
                
            }
            indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
            vertices.push_back({vec3(start_pos[0], start_pos[1], 0), vec3(1), vec3(1), vec4(0, 0, 0, 1)});
            
            if(propShowPoints)
            {
                indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(start_pos[0], start_pos[1], 0), vec3(1), vec3(1), vec4(0, 0, 0, 1)});
            }

        }
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = vectorField.getBBoxMin();
    const dvec2 bBoxMax = vectorField.getBBoxMax();

  
    // The start point should be inside the volume (set maximum to the upper
    // right corner)
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propSeedMode.get() == 0) {

        // Draw start point
        vec2 start_pos = propStartPoint.get();
        drawStreamline(mesh,vectorField,vertices,start_pos);
        // linspace(1, 3, 4); 
        vector<double> x_sp = linspace(-1, 1, 4);
        LogProcessorInfo(" X_sp.size() -> " << x_sp.size());
        for(int i = 0; i < x_sp.size(); i++)
        {
            LogProcessorInfo(" Linspace -> " << x_sp[i]);
        }


    } else {

        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
     // Extract the minimum and maximum value from the input data
   // const double minValue = grid.getMinValue();
    //const double maxValue = grid.getMaxValue();

        // Seed randomly Task 4.3 (a)

        
        double min_x = bBoxMin[0];
        double min_y = bBoxMin[1];
        double max_x =  bBoxMax[0];
        double max_y = bBoxMax[1];

        if(propUniformGrid)
        {
            // Uniform grid Task 4.3 (b)
            vec2 theGridPoints = propGridpoints.get();
            vector<double> X_axis = linspace(min_x, max_x, theGridPoints[0]);
            vector<double> Y_axis = linspace(min_y, max_y, theGridPoints[1]);

            for(int i = 0; i < X_axis.size(); i++)
            {
                for(int j = 0; j < Y_axis.size(); j++)
                {
                    // uniform grid 
                    dvec2 seedpoint = dvec2(X_axis[i],Y_axis[j]);
                    drawStreamline(mesh,vectorField,vertices,seedpoint);
                }
            }
        }

        if(propRandomGrid)
        {
            for(int n = 0; n < propNumberRamdomPoints; n++)
            {
               // dvec2 randomPoint = generateRandom(bBoxMin, bBoxMax);
               	//int seed = 1999;
                //srand(seed);
                double number_x = ((double)rand() / RAND_MAX) * (max_x - min_x) + min_x;
                double number_y = ((double)rand() / RAND_MAX) * (max_y - min_y) + min_y;
                dvec2 randomPoint = dvec2(number_x, number_y);

                dvec2 new_vec = vectorField.interpolate(randomPoint);
                LogProcessorInfo(" index -> " << n << " , vector (Rnd)-> " << randomPoint );
                LogProcessorInfo(" index -> " << n << " , vector (nVec)-> " << new_vec );

                drawStreamline(mesh,vectorField,vertices,randomPoint);

            }
        }

    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

}  // namespace inviwo
