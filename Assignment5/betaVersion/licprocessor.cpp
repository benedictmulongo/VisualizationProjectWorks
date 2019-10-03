/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>
#include <math.h>
#include <vector> 
#include <random>
#include <cstdlib> // required for srand(), rand().
using namespace std;
#define LINE_SQUARE_CLIP_MAX 100000

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
    , propLicMethods("licMethods", "Lic methods", 0, 0, 2, 1)
    , propKernelSize("kernelSize", "Kernel Size", 3, 2, 20, 1)
    , propContrastControl("contrastControl", "Contrast Control", 2, 1, 15, 1)
    , propContrast("contrast", "Contrast", 0, 0, 2, 1)
// TODO: Register additional properties
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
   // TODO: Register additional properties
    addProperty(propLicMethods);
    addProperty(propKernelSize);
    addProperty(propContrastControl);
    addProperty(propContrast);

}

double LICProcessor::velocity(const VectorField2& vectorField, const dvec2& pos)
{
    dvec2 new_vec = vectorField.interpolate(pos);
    double length = sqrt(new_vec[0]*new_vec[0] + new_vec[1]* new_vec[1] );
    return length;
}


vector<dvec2> LICProcessor::streamliners(VectorField2 vectorField, dvec2 start_pos, double propStepsize, bool propIntegration, const size2_t texDims_ ,  const size3_t vectorFieldDims_)
{
    vector<dvec2> curves;
    double propMaxArcLength = 30;
    int propMaxIntegrationSteps = 50;
    int Line_length = propKernelSize;

    curves.push_back(start_pos);

    // This is the maximun of arc length Task. 4.2 (e)
    float streamline_length = 0;

    for (int i = 0; i < Line_length; i++)
    {
        // Stop the integration after N stepsize Task 4.2 (d)
        if(i > propMaxIntegrationSteps) break; 

        double stepdirection = (propIntegration == true) ? propStepsize : -propStepsize;
        //start_pos = Integrator::RK4_lic(vectorField, start_pos, stepdirection, texDims_, vectorFieldDims_);
        start_pos = Integrator::RK4_stream(vectorField, start_pos, stepdirection);

        // Stop the integration when the velocity become slow Task 4.2 (h)
        if(velocity(vectorField, start_pos) < 1.0e-5) break;

        // Stop the integration at the zeros of the vector field Task 4.2 (g)
        if(start_pos[0]== 0 && start_pos[1] == 0) break;

        curves.push_back(start_pos);
    }

   /* for(int i = 0; streamline_length < propMaxArcLength; ++i, streamline_length += propStepsize)
    {

        // Stop the integration after N stepsize Task 4.2 (d)
        if(i > propMaxIntegrationSteps) break; 

        double stepdirection = (propIntegration == true) ? propStepsize : -propStepsize;
        //start_pos = Integrator::RK4_lic(vectorField, start_pos, stepdirection, texDims_, vectorFieldDims_);
        start_pos = Integrator::RK4_stream(vectorField, start_pos, stepdirection);

        // Stop the integration when the velocity become slow Task 4.2 (h)
        if(velocity(vectorField, start_pos) < 1.0e-5) break;

        // Stop the integration at the zeros of the vector field Task 4.2 (g)
        if(start_pos[0]== 0 && start_pos[1] == 0) break;

        curves.push_back(start_pos);
    }*/

    return curves;
}

double LICProcessor::kernelBox(const int t, const double mini, const double maxi )
{
    if (t < mini || t > maxi)
        return 0.0;
    else
        return 1./ (maxi - mini);
}

bool LICProcessor::inBounds(const dvec2& p, const size2_t texDims_)
{
    return (p[0] >= 0 && p[0] < texDims_.x && p[1] >= 0 && p[1] < texDims_.y);
}

double LICProcessor::LICPoint(const VectorField2& vField, const dvec2& posit, const RGBAImage text, const size2_t texDims_, const size3_t vectorFieldDims_)
{
	//map texture position to vector field position
	/*double pxWidth = (double)(vectorFieldDims_.x - 1) / (texDims_.x-1);
    double pxHeight = (double)( vectorFieldDims_.y - 1) / (texDims_.y - 1);
	double step = pxWidth > pxHeight ? pxWidth: pxHeight;*/

    double pxWidth = (double)(vField.getBBoxMax().x - vField.getBBoxMin().x) / texDims_.x;
    double pxHeight = (double)(vField.getBBoxMax().y - vField.getBBoxMin().y) / texDims_.y;
    double step = pxWidth > pxHeight ? pxWidth : pxHeight;

    dvec2 vf_bboxmin = vField.getBBoxMin();
    dvec2 scaleFactors = {
        (texDims_.x - 1) / (vField.getBBoxMax().x - vField.getBBoxMin().x),
        (texDims_.y - 1) / (vField.getBBoxMax().y - vField.getBBoxMin().y)};
    dvec2 scaleSlow = {
        (vField.getBBoxMax().x - vField.getBBoxMin().x) / (texDims_.x - 1),
        (vField.getBBoxMax().y - vField.getBBoxMin().y) / (texDims_.y - 1)};

    double mapX = ((double) posit[0]) * scaleSlow.x + vf_bboxmin.x;
    double mapY = ((double) posit[1]) * scaleSlow.y + vf_bboxmin.y;
    dvec2 pos = dvec2(mapX, mapY);

    //dvec2 pos = dvec2(posit[0] * pxWidth, posit[1] * pxHeight );

    vector<dvec2> curv_forward =  streamliners(vField, pos, step, true,texDims_, vectorFieldDims_); 
    vector<dvec2> curv_backward = streamliners(vField, pos, step, false,texDims_, vectorFieldDims_); 

    double valuer = 0;
	for (int t = curv_backward.size()-1; t >= 1; t--) {
        double textVal = 0 ;
		auto pt = curv_backward[t];
        int mapPosX = (curv_backward[t].x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
               (vField.getBBoxMax().x - vField.getBBoxMin().x);
        int mapPosY = (curv_backward[t].y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
               (vField.getBBoxMax().y - vField.getBBoxMin().y);

		//map vector field position to texture position
		//int mapPosX = int(pt.x * (double)(texDims_.x - 1)/ (vectorFieldDims_.x - 1));
		//int mapPosY = int(pt.y * (double)(texDims_.y - 1) / (vectorFieldDims_.y - 1));
        
		//auto color = colLayer->getAsDVec4(size2_t(mapPosX, mapPosY));
        
        if( inBounds(dvec2(mapPosX, mapPosY), texDims_) )
        {
            textVal = text.readPixelGrayScale(size2_t(mapPosX, mapPosY));
        }
		auto k = kernelBox(t, 0, curv_backward.size() + curv_forward.size());
		valuer += textVal * k;
	}

	for (int t = 0; t < curv_forward.size(); t++) {
        double textVal = 0 ;
		auto pt = curv_forward[t];

        int mapPosX = (curv_forward[t].x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
               (vField.getBBoxMax().x - vField.getBBoxMin().x);
        int mapPosY = (curv_forward[t].y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
               (vField.getBBoxMax().y - vField.getBBoxMin().y);

		//map vector field position to texture position
		//int mapPosX = int(pt.x * (double)(texDims_.x - 1)/ (vectorFieldDims_.x - 1));
		//int mapPosY = int(pt.y * (double)(texDims_.y - 1) / (vectorFieldDims_.y - 1));
		//auto color = colLayer->getAsDVec4(size2_t(mapPosX, mapPosY));
        if( inBounds(dvec2(mapPosX, mapPosY), texDims_) )
        {
            textVal = text.readPixelGrayScale(size2_t(mapPosX, mapPosY));
        }
		auto k = kernelBox(t, 0, curv_backward.size() + curv_forward.size());
		valuer += textVal * k;
	}

    //valuer = valuer / ((curv_forward.size() +  curv_backward.size() + 1)*3) ;
    //valuer = valuer / (curv_forward.size() +  curv_backward.size() + 1) ;

    return valuer; 
}

dvec2 LICProcessor::getIndices(const int row, const int col, const int index)
{
    int count = 0;
    dvec2 indexes = dvec2(0,0);
    for(int y = 0; y < col; y++)
    {
        for(int x = 0; x < row; x++)
        {
            if(count == index)
            {
                indexes[0] = y;
                indexes[1] = x;
                break; 
            }
            count = count + 1;
        } 
    }
    return indexes;
}

vector<vector<double>> LICProcessor::LICAlgo(const VectorField2& vField, const RGBAImage text, const size2_t texDims_, const size3_t vectorFieldDims_)
{
    vector<vector<double>> mat(texDims_.x, vector<double>(texDims_.y) ); 

	for (int j = 0; j < texDims_.y; j++) 
    {
		for (int i = 0; i < texDims_.x; i++)
        {
			dvec2 position = dvec2(i,j);
            mat[i][j] = LICPoint(vField, position, text, texDims_, vectorFieldDims_);
			//mat[i][j] = LICPoint(tr,texDims,vr,dims, position);	
		}
	}

	return mat;
}

vector<vector<double>> LICProcessor::FastLic(const VectorField2& vField, const RGBAImage text, const size2_t texDims_, const size3_t vectorFieldDims_)
{
    vector<vector<double>> mat(texDims_.x, vector<double>(texDims_.y) ); 
    vector<vector<double>> numHits(texDims_.x, vector<double>(texDims_.y, 0.0) ); 
    vector<vector<double>> accumInt(texDims_.x, vector<double>(texDims_.y, 0.0) ); 

    int L = 20;
    double M = 10;
    double minNumHits = 1; 
    double m = 0;

	for (int j = 0; j < texDims_.y; j++) 
    {
		for (int i = 0; i < texDims_.x; i++)
        {
            if(numHits[i][j] < minNumHits) 
            {
                dvec2 posit = dvec2(i + 0.5,j + 0.5);

                double pxWidth = (double)(vField.getBBoxMax().x - vField.getBBoxMin().x) / texDims_.x;
                double pxHeight = (double)(vField.getBBoxMax().y - vField.getBBoxMin().y) / texDims_.y;
                double step = pxWidth > pxHeight ? pxWidth : pxHeight;

                dvec2 vf_bboxmin = vField.getBBoxMin();
                dvec2 scaleFactors = {
                    (texDims_.x - 1) / (vField.getBBoxMax().x - vField.getBBoxMin().x),
                    (texDims_.y - 1) / (vField.getBBoxMax().y - vField.getBBoxMin().y)};
                dvec2 scaleSlow = {
                    (vField.getBBoxMax().x - vField.getBBoxMin().x) / (texDims_.x - 1),
                    (vField.getBBoxMax().y - vField.getBBoxMin().y) / (texDims_.y - 1)};

                double mapX = ((double) posit[0]) * scaleSlow.x + vf_bboxmin.x;
                double mapY = ((double) posit[1]) * scaleSlow.y + vf_bboxmin.y;
                dvec2 pos = dvec2(mapX, mapY);

                vector<dvec2> curv_forward =  streamliners(vField, pos, step, true,texDims_, vectorFieldDims_); 
                vector<dvec2> curv_backward = streamliners(vField, pos, step, false,texDims_, vectorFieldDims_); 

                double valuer = 0;
                for (int t = curv_backward.size()-1; t >= 1; t--) {
                    double textVal = 0 ;
                    auto pt = curv_backward[t];
                    int mapPosX = (curv_backward[t].x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
                        (vField.getBBoxMax().x - vField.getBBoxMin().x);
                    int mapPosY = (curv_backward[t].y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
                        (vField.getBBoxMax().y - vField.getBBoxMin().y);
                    
                    if( inBounds(dvec2(mapPosX, mapPosY), texDims_) )
                    {
                        textVal = text.readPixelGrayScale(size2_t(mapPosX, mapPosY));
                    }
                    auto k = kernelBox(t, 0, curv_backward.size() + curv_forward.size());
                    valuer += textVal * k;
                }

                for (int t = 0; t < curv_forward.size(); t++) {
                    double textVal = 0 ;
                    auto pt = curv_forward[t];

                    int mapPosX = (curv_forward[t].x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
                        (vField.getBBoxMax().x - vField.getBBoxMin().x);
                    int mapPosY = (curv_forward[t].y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
                        (vField.getBBoxMax().y - vField.getBBoxMin().y);

                    if( inBounds(dvec2(mapPosX, mapPosY), texDims_) )
                    {
                        textVal = text.readPixelGrayScale(size2_t(mapPosX, mapPosY));
                    }
                    auto k = kernelBox(t, 0, curv_backward.size() + curv_forward.size());
                    valuer += textVal * k;
                }

                vector<double> Intensity; 
                Intensity.assign((curv_backward.size() + curv_forward.size()), 0);
                vector<double> Intensity_b; 
                Intensity_b.assign((curv_backward.size() + curv_forward.size()), 0);

                mat[i][j] = valuer;
                m = 1;
                Intensity[0] = valuer;
                Intensity_b[0] = valuer;
                while(m < ((curv_backward.size() + curv_forward.size() ) / 2 ) )
                {
                    if (m < curv_forward.size() && (curv_forward.size()-m > 0) ) 
                    {
                        dvec2 forward1 = curv_forward[m]; // x_i + (Nf + 1)
                        dvec2 forward2 = curv_forward[curv_forward.size()-m]; // 

                        forward1[0] = (forward1.x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
                            (vField.getBBoxMax().x - vField.getBBoxMin().x);
                        forward1[1] = (forward1.y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
                            (vField.getBBoxMax().y - vField.getBBoxMin().y);

                        forward2[0] = (forward2.x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
                            (vField.getBBoxMax().x - vField.getBBoxMin().x);
                        forward2[1] = (forward2.y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
                            (vField.getBBoxMax().y - vField.getBBoxMin().y);

                        auto k = kernelBox(m, 0, curv_backward.size() + curv_forward.size());

                        // Put at point forward[m]
                        Intensity[m] = Intensity[m-1] - (text.readPixelGrayScale(size2_t(forward2[0], forward2[1])) * k ) + (text.readPixelGrayScale(size2_t(forward1[0], forward1[1])) * k ) ;
                        mat[(int)forward1[0]][(int)forward1[1]] = Intensity[m];
                        numHits[(int)forward1[0]][(int)forward1[1]]  += 1;

                    }

                    if (m < curv_backward.size() && (curv_backward.size()-m > 0) ) 
                    {
                        dvec2 backward1 = curv_backward[m]; // x_i + (Nf + 1)
                        dvec2 backward2 = curv_backward[curv_backward.size()-m]; // 

                        backward1[0] = (backward1.x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
                            (vField.getBBoxMax().x - vField.getBBoxMin().x);
                        backward1[1] = (backward1.y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
                            (vField.getBBoxMax().y - vField.getBBoxMin().y);

                        backward2[0] = (backward2.x - vField.getBBoxMin().x) * (double)(texDims_.x - 1) /
                            (vField.getBBoxMax().x - vField.getBBoxMin().x);
                        backward2[1] = (backward2.y - vField.getBBoxMin().y) * (double)(texDims_.y - 1) /
                            (vField.getBBoxMax().y - vField.getBBoxMin().y);

                        auto k = kernelBox(m, 0, curv_backward.size() + curv_forward.size());

                        // Put at point forward[m]
                        Intensity_b[m] = Intensity_b[m-1] - (text.readPixelGrayScale(size2_t(backward2[0], backward2[1])) * k ) + (text.readPixelGrayScale(size2_t(backward1[0], backward1[1])) * k ) ;
                        mat[(int)backward1[0]][(int)backward1[1]] = Intensity[m];
                        numHits[(int)backward1[0]][(int)backward1[1]] += 1;

                    }

                    m = m + 1 ; 
                }
            }

		}
	}

	return mat;

}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    //const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    // TODO: Implement LIC and FastLIC
    // This code instead just creates a black image
    //LIC(vectorField, licImage, texture, texDims_);
    //licTexture = LICAlgo(vectorField, texture, texDims_, vectorFieldDims_);
    //licTexture = FastLic(vectorField, texture, texDims_, vectorFieldDims_);

    if(propLicMethods == 0)
    {
        for (auto j = 0; j < texDims_.y; j++)
        {
            for (auto i = 0; i < texDims_.x; i++)
            {
                double val = texture.readPixelGrayScale(size2_t(i, j));
                licImage.setPixelGrayScale(size2_t(i,j), val);
            }
        }
    }
    else if(propLicMethods == 1)
    {
        // FastLic 
        licTexture = FastLic(vectorField, texture, texDims_, vectorFieldDims_);

        if(propContrast == 0)
        {
            for(auto j = 0; j < texDims_.y; j++)
            {
                for(auto i = 0; i < texDims_.x; i++)
                {
                    double val = licTexture[i][j];
                    licImage.setPixelGrayScale(size2_t(i,j), val);
                }
            }
        }
        else
        {
            for(auto j = 0; j < texDims_.y; j++)
            {
                for(auto i = 0; i < texDims_.x; i++)
                {
                    
                    double val1 = (inBounds(dvec2(i + 1,j - 1), texDims_)) ? licTexture[i + 1][j - 1] : 0;
                    double val2 = (inBounds(dvec2(i + 1,j), texDims_)) ? licTexture[i + 1][j] : 0;
                    double val3 = (inBounds(dvec2(i + 1,j + 1), texDims_)) ? licTexture[i + 1][j + 1] : 0;
                    double val4 = (inBounds(dvec2(i,j + 1), texDims_)) ? licTexture[i][j + 1] : 0;
                    double val5 = (inBounds(dvec2(i - 1,j + 1), texDims_)) ? licTexture[i - 1][j + 1] : 0;
                    double val6 = (inBounds(dvec2(i - 1,j), texDims_)) ? licTexture[i - 1][j] : 0;
                    double val7 = (inBounds(dvec2(i - 1,j - 1), texDims_)) ? licTexture[i - 1][j - 1] : 0;
                    double val8 = (inBounds(dvec2(i,j - 1), texDims_)) ? licTexture[i][j - 1] : 0;

                    double mean = (val1 + val2 + val3 + val4 + val5 + val6 + val7 + val8) / 8 ;
                    double std1 = (val1 - mean) * (val1 - mean);
                    double std2 = (val2 - mean) * (val2 - mean);
                    double std3 = (val3 - mean) * (val3 - mean);
                    double std4 = (val4 - mean) * (val4 - mean);
                    double std5 = (val5 - mean) * (val5 - mean);
                    double std6 = (val6 - mean) * (val6 - mean);
                    double std7 = (val7 - mean) * (val7 - mean);
                    double std8 = (val8 - mean) * (val8 - mean);

                    double std = ( std1 + std2 + std3 + std4 + std5 + std6 + std7 + std8) / 8 ;

                    double f_ij = mean + ( propContrastControl / std )* (licTexture[i][j] - mean) ;
    


                    licImage.setPixelGrayScale(size2_t(i,j), f_ij);
                }
            }
        }
        

    }
    else
    {
        // Slow Lic
        //licTexture = LICAlgo(vectorField, texture, texDims_, vectorFieldDims_);
        for (auto j = 0; j < texDims_.y; j++)
        {
            for (auto i = 0; i < texDims_.x; i++)
            {
                //double val = licTexture[i][j];
                //licImage.setPixelGrayScale(size2_t(i,j), val);

                dvec2 position = dvec2(i,j);
                double val = LICPoint(vectorField, position, texture, texDims_, vectorFieldDims_);
                licImage.setPixelGrayScale(size2_t(i,j), val);

            }
        }
    }
    


    licOut_.setData(outImage);
}

}  // namespace inviwo
