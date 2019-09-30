/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>
#include <cmath>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",  // Class identifier
    "Marching Squares",            // Display name
    "KTH Lab",                     // Category
    CodeState::Experimental,       // Code state
    Tags::None,                    // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshOut("meshOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData) {
    // Register ports
    addPort(inData);
    addPort(meshOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("midpoint", "Mid Point", 0);
    propDeciderType.addOption("asymptotic", "Asymptotic", 1);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propNumContours, propIsoTransferFunc);

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            // util::hide(propIsoValue);
            // util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            util::hide(propIsoValue, propIsoColor);
            util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    double minValue = grid.getMinValue();
    double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    // LogProcessorInfo("This scalar field contains values between " << minValue << " and " <<
    // maxValue << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    // LogProcessorInfo("The value at (0,0) is: " << valueAt00 << ".");

    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        // TODO: Add grid lines of the given color

        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
        auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        int x_vertices = nVertPerDim[0];
        int y_vertices = nVertPerDim[1];
        float x_cellSize = cellSize[0];
        float y_cellSize = cellSize[1];

        for (size_t j = 0; j < y_vertices; j++) {
            vec2 v1 = vec2(0, j * y_cellSize);
            vec2 v2 = vec2(x_cellSize * (x_vertices - 1), j * y_cellSize);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), vertices);
        }
        for (size_t i = 0; i < x_vertices; i++) {
            vec2 v1 = vec2(i * x_cellSize, 0);
            vec2 v2 = vec2(i * x_cellSize, y_cellSize * (y_vertices - 1));
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), vertices);
        }
    }  // LogProcessorInfo("CELLSIZE: " << cellSize);

    // Iso contours

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    int x_vertices = nVertPerDim[0];
    int y_vertices = nVertPerDim[1];

    /*ScalarField2 smoothedField = ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin);
    int sigma = 2;
    
    const double PI = 3.141592653589793238463;
    for (int i = 0; i < x_vertices; i++) {
        for (int j = 0; j < y_vertices; j++) {
                        float gFilter = (1/(2*PI*pow(sigma,
    2)))*(exp(-((pow(i,2)+pow(j,2))/(2*pow(sigma, 2))))); if (i == 0 && j == 0) minValue = gFilter;
            smoothedField.setValueAtVertex({i, j}, gFilter);
            if (gFilter < minValue) minValue = gFilter;
            //LogProcessorInfo("GFilter:" << gFilter);
                }
        }
    for (int i=0;i<x_vertices;i++)
    //minValue = smoothedField.getMinValue();
    maxValue = smoothedField.getMaxValue();
    //LogProcessorInfo("Min:" << minValue << " Max:" << maxValue);
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);*/

    auto indexBufferIsoLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    if (propMultiple.get() == 0) {
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
        // and color it with the specified color (propIsoColor)
        // do something
        // LogProcessorInfo("CELLSIZE: " << cellSize);
        for (size_t i = 0; i < nVertPerDim[0] - 1; i++) {
            for (size_t j = 0; j < nVertPerDim[1] - 1; j++) {
                float value_A = grid.getValueAtVertex({i, j});
                float value_B = grid.getValueAtVertex({i, j + 1});
                float value_C = grid.getValueAtVertex({i + 1, j + 1});
                float value_D = grid.getValueAtVertex({i + 1, j});

                std::array<bool, 4> corners;
                std::array<bool, 4> bipolar;
                std::array<vec2, 4> edges;

                bipolar[0] = (propIsoValue >= value_A && propIsoValue <= value_B) ||
                             (propIsoValue >= value_B && propIsoValue <= value_A);
                corners[0] = value_A <= propIsoValue;
                bipolar[1] = (propIsoValue >= value_B && propIsoValue <= value_C) ||
                             (propIsoValue >= value_C && propIsoValue <= value_B);
                corners[1] = value_B <= propIsoValue;
                bipolar[2] = (propIsoValue >= value_C && propIsoValue <= value_D) ||
                             (propIsoValue >= value_D && propIsoValue <= value_C);
                corners[2] = value_C <= propIsoValue;
                bipolar[3] = (propIsoValue >= value_D && propIsoValue <= value_A) ||
                             (propIsoValue >= value_A && propIsoValue <= value_D);
                corners[3] = value_D <= propIsoValue;

                // LogProcessorInfo("A: " << value_A << " B: " << value_B << " C: " << value_C << "
                // D: " << value_D);

                // calculating cuttingPoints if necessary
                for (size_t k = 0; k < corners.size(); k++) {
                    float high_value;
                    float low_value;
                    float cuttingPoint;

                    if (bipolar[k]) {
                        switch (k) {
                            case 0:
                                high_value = value_B;
                                low_value = value_A;
                                break;
                            case 1:
                                high_value = value_C;
                                low_value = value_B;
                                break;
                            case 2:
                                high_value = value_C;
                                low_value = value_D;
                                break;
                            case 3:
                                high_value = value_D;
                                low_value = value_A;
                                break;
                            default:
                                break;
                        }

                        cuttingPoint = (propIsoValue - low_value) / (high_value - low_value);

                        switch (k) {
                            case 0:
                                edges[k] = vec2(i * cellSize[0], (j + cuttingPoint) * cellSize[1]);
                                break;
                            case 1:
                                edges[k] =
                                    vec2((i + cuttingPoint) * cellSize[0], (j + 1) * cellSize[1]);
                                break;
                            case 2:
                                edges[k] =
                                    vec2((i + 1) * cellSize[0], (j + cuttingPoint) * cellSize[1]);
                                break;
                            case 3:
                                edges[k] = vec2((i + cuttingPoint) * cellSize[0], j * cellSize[1]);
                                break;
                            default:
                                break;
                        }
                    }
                }

                /*LogProcessorInfo("CORNERS: " << corners[0] << " " << corners[1] << " " <<
                corners[2]
                                             << " " << corners[3]);
                LogProcessorInfo("EDGES: " << edges[0] << " " << edges[1] << " " << edges[2] << " "
                                           << edges[3]);
                LogProcessorInfo("BIPOLAR: " << bipolar[0] << " " << bipolar[1] << " " << bipolar[2]
                                             << " " << bipolar[3]); */

                int bipolarCounter = 0;
                int cornerCounter = 0;
                for (size_t k = 0; k < 4; k++) {
                    if (bipolar[k]) {
                        bipolarCounter++;
                    }
                    if (corners[k]) {
                        cornerCounter++;
                    }
                }

                vec2 edge1;
                vec2 edge2;
                vec2 edge3;
                vec2 edge4;
                std::vector<vec2> two_edges = std::vector<vec2>();
                switch (bipolarCounter) {
                    case 4: {

                        // midpoint decider
                        if ((1 / 4) * (value_A + value_B + value_C + value_D) <= propIsoValue) {
                            if (corners[0]) {
                                edge1 = edges[0];
                                edge2 = edges[3];
                                edge3 = edges[1];
                                edge4 = edges[2];
                            } else {
                                edge1 = edges[0];
                                edge2 = edges[1];
                                edge3 = edges[2];
                                edge4 = edges[3];
                            }
                        } else {
                            if (corners[0]) {
                                edge1 = edges[0];
                                edge2 = edges[3];
                                edge3 = edges[1];
                                edge4 = edges[2];
                            } else {
                                edge1 = edges[0];
                                edge2 = edges[1];
                                edge3 = edges[2];
                                edge4 = edges[3];
                            }
                        }

                        drawLineSegment(edge1, edge2, propIsoColor, indexBufferIsoLines.get(),
                                        vertices);
                        drawLineSegment(edge3, edge4, propIsoColor, indexBufferIsoLines.get(),
                                        vertices);
                        break;
                    }
                    case 2: {
                        for (size_t k = 0; k < bipolar.size(); k++) {
                            if (bipolar[k]) {
                                two_edges.push_back(edges[k]);
                            }
                        }
                        drawLineSegment(two_edges.at(0), two_edges.at(1), propIsoColor,
                                        indexBufferIsoLines.get(), vertices);
                    }
                    default:
                        break;
                }
            }
        }

    } else {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value
        LogProcessorInfo("Num of contours:" << propNumContours);
        LogProcessorInfo("Min:" << minValue << " Max:" << maxValue);
        float stepSize = (maxValue - minValue) / propNumContours;
        for (float iso = minValue; iso <= maxValue; iso += stepSize) {
            LogProcessorInfo("Iso:" << iso);
            for (size_t i = 0; i < nVertPerDim[0] - 1; i++) {
                for (size_t j = 0; j < nVertPerDim[1] - 1; j++) {
                    float value_A = grid.getValueAtVertex({i, j});
                    float value_B = grid.getValueAtVertex({i, j + 1});
                    float value_C = grid.getValueAtVertex({i + 1, j + 1});
                    float value_D = grid.getValueAtVertex({i + 1, j});

                    std::array<bool, 4> corners;
                    std::array<bool, 4> bipolar;
                    std::array<vec2, 4> edges;

                    bipolar[0] =
                        (iso >= value_A && iso <= value_B) || (iso >= value_B && iso <= value_A);
                    corners[0] = value_A <= iso;
                    bipolar[1] =
                        (iso >= value_B && iso <= value_C) || (iso >= value_C && iso <= value_B);
                    corners[1] = value_B <= iso;
                    bipolar[2] =
                        (iso >= value_C && iso <= value_D) || (iso >= value_D && iso <= value_C);
                    corners[2] = value_C <= iso;
                    bipolar[3] =
                        (iso >= value_D && iso <= value_A) || (iso >= value_A && iso <= value_D);
                    corners[3] = value_D <= iso;

                    // LogProcessorInfo("A: " << value_A << " B: " << value_B << " C: " << value_C
                    // << " D: " << value_D);

                    // calculating cuttingPoints if necessary
                    for (size_t k = 0; k < corners.size(); k++) {
                        float high_value;
                        float low_value;
                        float cuttingPoint;

                        if (bipolar[k]) {
                            switch (k) {
                                case 0:
                                    high_value = value_B;
                                    low_value = value_A;
                                    break;
                                case 1:
                                    high_value = value_C;
                                    low_value = value_B;
                                    break;
                                case 2:
                                    high_value = value_C;
                                    low_value = value_D;
                                    break;
                                case 3:
                                    high_value = value_D;
                                    low_value = value_A;
                                    break;
                                default:
                                    break;
                            }

                            cuttingPoint = (iso - low_value) / (high_value - low_value);

                            switch (k) {
                                case 0:
                                    edges[k] =
                                        vec2(i * cellSize[0], (j + cuttingPoint) * cellSize[1]);
                                    break;
                                case 1:
                                    edges[k] = vec2((i + cuttingPoint) * cellSize[0],
                                                    (j + 1) * cellSize[1]);
                                    break;
                                case 2:
                                    edges[k] = vec2((i + 1) * cellSize[0],
                                                    (j + cuttingPoint) * cellSize[1]);
                                    break;
                                case 3:
                                    edges[k] =
                                        vec2((i + cuttingPoint) * cellSize[0], j * cellSize[1]);
                                    break;
                                default:
                                    break;
                            }
                        }
                    }

                    /*LogProcessorInfo("CORNERS: " << corners[0] << " " << corners[1] << " " <<
                    corners[2]
                                                                             << " " << corners[3]);
                    LogProcessorInfo("EDGES: " << edges[0] << " " << edges[1] << " " << edges[2] <<
                    " "
                                                                       << edges[3]);
                    LogProcessorInfo("BIPOLAR: " << bipolar[0] << " " << bipolar[1] << " " <<
                    bipolar[2]
                                                                             << " " << bipolar[3]);
                  */

                    int bipolarCounter = 0;
                    int cornerCounter = 0;
                    for (size_t k = 0; k < 4; k++) {
                        if (bipolar[k]) {
                            bipolarCounter++;
                        }
                        if (corners[k]) {
                            cornerCounter++;
                        }
                    }

                    vec2 edge1;
                    vec2 edge2;
                    vec2 edge3;
                    vec2 edge4;
                    float normalizedIso = float((iso - minValue) / (maxValue - minValue));
                    // LogProcessorInfo("Normalized iso:" << normalizedIso);
                    vec4 color = propIsoTransferFunc.get().sample(normalizedIso);
                    // LogProcessorInfo("Color red:" << color.r << " green:" << color.g << " blue:"
                    // << color.b);
                    std::vector<vec2> two_edges = std::vector<vec2>();
                    switch (bipolarCounter) {
                        case 4: {

                            // midpoint decider
                            if ((1 / 4) * (value_A + value_B + value_C + value_D) <= iso) {
                                if (corners[0]) {
                                    edge1 = edges[0];
                                    edge2 = edges[3];
                                    edge3 = edges[1];
                                    edge4 = edges[2];
                                } else {
                                    edge1 = edges[0];
                                    edge2 = edges[1];
                                    edge3 = edges[2];
                                    edge4 = edges[3];
                                }
                            } else {
                                if (corners[0]) {
                                    edge1 = edges[0];
                                    edge2 = edges[3];
                                    edge3 = edges[1];
                                    edge4 = edges[2];
                                } else {
                                    edge1 = edges[0];
                                    edge2 = edges[1];
                                    edge3 = edges[2];
                                    edge4 = edges[3];
                                }
                            }

                            drawLineSegment(edge1, edge2, color, indexBufferIsoLines.get(),
                                            vertices);
                            drawLineSegment(edge3, edge4, color, indexBufferIsoLines.get(),
                                            vertices);
                            break;
                        }
                        case 2: {
                            for (size_t k = 0; k < bipolar.size(); k++) {
                                if (bipolar[k]) {
                                    two_edges.push_back(edges[k]);
                                }
                            }
                            drawLineSegment(two_edges.at(0), two_edges.at(1), color,
                                            indexBufferIsoLines.get(), vertices);
                        }
                        default:
                            break;
                    }
                }
            }
        }

        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo