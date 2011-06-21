/*!
 * @author  Xinghua Lou <xinghua.lou@iwr.uni-heidelberg.de>
 *
 * @section LICENSE
 * 
 * BOT. Copyright (c) 2010 by Xinghua Lou.
 *
 * This software was developed by Xinghua Lou.
 * Enquiries shall be directed to: xinghua.lou@iwr.uni-heidelberg.de.
 *
 * All advertising materials mentioning features or use of this software must
 * display the following acknowledgement: ``This product includes the BOT
 * library developed by Xinghua Lou. Please direct enquiries concerning BOT to 
 * xinghua.lou@iwr.uni-heidelberg.de.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, 
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - All advertising materials mentioning features or use of this software must 
 *   display the following acknowledgement: ``This product includes the BOT
 *   library developed by Xinghua Lou. Please direct enquiries concerning BOT to 
 *   xinghua.lou@iwr.uni-heidelberg.de.
 * - The names of the authors must not be used to endorse or promote products 
 *   derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
 * EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "ObjectFeatureFactory.hxx"
#include "SingletsGenerator.hxx"
#include "MultipletsGenerator.hxx"
#include "InputOutput.hxx"
#include "vigra/hdf5impex.hxx"
#include <vector>

using namespace bot;

int main()
{
    // test on real data - 1
    {
        std::cout << "**** Test on real data (time = 1) ****" << std::endl;
        // load from hdf5 file
        vigra::HDF5ImportInfo info_img("../data/test-image-seg-1.h5", "/raw/volume");
        vigra::MultiArrayShape<2 >::type shape_img(info_img.shape()[0], info_img.shape()[1]);
        vigra::MultiArray<2, unsigned short > img_(shape_img, static_cast<unsigned short >(0));
        vigra::readHDF5(info_img, img_);
        Matrix2D img(img_);

        vigra::HDF5ImportInfo info_seg("../data/test-image-seg-1.h5", "/segmentation/volume");
        vigra::MultiArrayShape<2 >::type shape_seg(info_seg.shape()[0], info_seg.shape()[1]);
        vigra::MultiArray<2, unsigned short > seg_(shape_seg, static_cast<unsigned short >(0));
        vigra::readHDF5(info_seg, seg_);
        Matrix2D seg(seg_);

        std::cout << "Read in image and segementation: " << std::endl <<
            "\timage shape = " << img.shape() << std::endl << 
            "\tsegmentation shape = " << seg.shape() << std::endl;

        // create the object list
        SingletsGenerator singlets_gen;
        Singlets singlets = singlets_gen(img, seg);
        std::cout << "Generate object list: " << singlets.size() << " objects" << std::endl;
        std::cout << singlets << std::endl;

        // create the party list
        MultipletsGenerator multiplets_gen(4, 40);
        Multiplets multiplets = multiplets_gen(img, seg, singlets);
        std::cout << "Generate party list: " << multiplets.size() << " parties" << std::endl;
        std::cout << multiplets << std::endl;
    }

    // test on real data - 2
    {
        std::cout << "**** Test on real data (time = 99) ****" << std::endl;
        // load from hdf5 file
        vigra::HDF5ImportInfo info_img("../data/test-image-seg-99.h5", "/raw/volume");
        vigra::MultiArrayShape<2 >::type shape_img(info_img.shape()[0], info_img.shape()[1]);
        vigra::MultiArray<2, unsigned short > img_(shape_img, static_cast<unsigned short >(0));
        vigra::readHDF5(info_img, img_);
        Matrix2D img(img_);

        vigra::HDF5ImportInfo info_seg("../data/test-image-seg-99.h5", "/segmentation/volume");
        vigra::MultiArrayShape<2 >::type shape_seg(info_seg.shape()[0], info_seg.shape()[1]);
        vigra::MultiArray<2, unsigned short > seg_(shape_seg, static_cast<unsigned short >(0));
        vigra::readHDF5(info_seg, seg_);
        Matrix2D seg(seg_);

        std::cout << "Read in image and segementation: " << std::endl <<
            "\timage shape = " << img.shape() << std::endl << 
            "\tsegmentation shape = " << seg.shape() << std::endl;

        // create the object list
        SingletsGenerator singlets_gen;
        Singlets singlets = singlets_gen(img, seg);
        std::cout << "Generate object list: " << singlets.size() << " objects" << std::endl;
        std::cout << singlets << std::endl;

        // create the party list
        MultipletsGenerator multiplets_gen(4, 40);
        Multiplets multiplets = multiplets_gen(img, seg, singlets);
        std::cout << "Generate party list: " << multiplets.size() << " parties" << std::endl;
        std::cout << multiplets << std::endl;
        
        // preapre the context
        Context context(img);
        std::cout << "Context::bounding_box = " << context.bounding_box() << std::endl;

        // dynamically create feature extractors
        std::string list[] = {
            "ObjectFeatureVolume", "ObjectFeatureSpatialDensityVariation", 
            "ObjectFeaturePosition", "ObjectFeatureBoundingBox", 
            "ObjectFeaturePrincipalComponents", "ObjectFeatureEccentricity", 
            "ObjectFeatureIntensityMean", "ObjectFeatureIntensityDeviation", 
            "ObjectFeatureIntensitySum", "ObjectFeatureIntensityHistogram",
            "ObjectFeatureBorderDistance", "ObjectFeatureBorderOverlap", 
            "ObjectFeatureShapeCompactness", "ObjectFeatureMassEvenness", 
            "ObjectFeatureVolumeEvenness", 
            "ObjectFeatureNull"
        };
        std::vector<ObjectFeatureFactory* > features;
        for (int i = 0; i < sizeof(list) / sizeof(std::string *); i ++) {
            ObjectFeatureFactory* fea = ObjectFeatureFactory::make(list[i]);
            if (!fea) 
                std::cerr << "Warning: skipping unknown feature class name " << list[i] << std::endl;
            else 
                features.push_back(fea);
        }

        {
            for (int i = 0; i < features.size(); i++) {
                std::cout << "**** Extracting " << list[i] << " ****" << std::endl;
                Matrix2D feature_mat;
                features[i]->extract(feature_mat, singlets.front(), context);

                std::cout << "Feature " << i << " for object " << singlets.front().id() << std::endl;
                std::cout << "\t\tShape: " << feature_mat.shape() << std::endl;
                std::cout << "\t\tValue: " << feature_mat << std::endl;
            }
        }

        std::cout << "-----------------------------" << std::endl;

        {
            for (int i = 0; i < features.size(); i++) {
                std::cout << "**** Extracting " << list[i] << " ****" << std::endl;
                Matrix2D feature_mat;
                features[i]->extract(feature_mat, multiplets.front(), context);

                std::cout << "Feature " << i << " for party " << multiplets.front().id() << std::endl;
                std::cout << "\t\tShape: " << feature_mat.shape() << std::endl;
                std::cout << "\t\tValue: " << feature_mat << std::endl;
            }
        }

        for (int i = 0; i < features.size(); i++)
            delete features[i];
    }

    return 0;
}
