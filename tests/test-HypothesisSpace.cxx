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

#include "InputOutput.hxx"
#include "HypothesisSpace.hxx"
#include "ObjectFeatureExtractor.hxx"
#include "AverageObject.hxx"
#include "SingletsGenerator.hxx"
#include "MultipletsGenerator.hxx"
#include "HDF5ReaderWriter.hxx"

using namespace bot;

void load_hdf5(const std::string& filename, 
               const std::string& image_path, 
               const std::string& segmentation_path, 
               std::vector<Matrix2D >& images,
               std::vector<Matrix2D >& segmentations)
{
    // load images
    vigra::HDF5ImportInfo info_img(filename.c_str(), image_path.c_str());
    vigra::MultiArrayShape<3 >::type shape_img(info_img.shape()[0], info_img.shape()[1], info_img.shape()[2]);
    vigra::MultiArray<3, uint16 > img_(shape_img, static_cast<unsigned short >(0));
    vigra::readHDF5(info_img, img_);
    for (int32 ind = 0; ind < img_.shape(2); ind ++) {
        Matrix2D img(img_.bindOuter(ind));
        images.push_back(img);
    }

    // load segmentations
    vigra::HDF5ImportInfo info_seg(filename.c_str(), segmentation_path.c_str());
    vigra::MultiArrayShape<3 >::type shape_seg(info_seg.shape()[0], info_seg.shape()[1], info_seg.shape()[2]);
    vigra::MultiArray<3, unsigned short > seg_(shape_seg, static_cast<unsigned short >(0));
    vigra::readHDF5(info_seg, seg_);
    for (int32 ind = 0; ind < seg_.shape(2); ind ++) {
        Matrix2D seg(seg_.bindOuter(ind));
        segmentations.push_back(seg);
    }
};

template<class T >
void print_features(const T& object)
{
    for (int32 ind = 0; ind < object.names().size(); ind ++) {
        std::cout << object.names()[ind] << ": " << object.features()[ind] << std::endl;
    }
};

int main()
{
    // load the image sequence
    std::string filename("../data/dcelliq-sequence-training.h5");
    std::vector<Matrix2D > images, segmentations;
    HDF5ReaderWriter::load(filename.c_str(), images, segmentations);

    // get the context
    Context context(images);
    std::cout << "****Context****" << std::endl << context << std::endl << std::endl;

    // load the configuration
    HypothesisSpace hypoSpace("../data/event-configuration-cell.ini");
    EventConfiguration conf = hypoSpace.configuration();

    // create singlets/muliplets and extract object features
    std::cout << "****Extracting singlets and multiplets" << std::endl;
    SingletsSequence singletsSequence;
    SingletsSequence avgSingletSequence;
    MultipletsSequence multipletsSequence;
    SingletsGenerator singletsGenerator;
    MultipletsGenerator multipletsGenerator(conf.k(), conf.d_max());
    ObjectFeatureExtractor objectFeatureExtractor(conf.get_feature_names(), context);
    for (int32 indT = 0; indT < images.size(); indT ++) {
        // generate singlets and multiplets
        Singlets singlets = singletsGenerator(images[indT], segmentations[indT]);
        Multiplets multiplets = multipletsGenerator(images[indT], segmentations[indT], singlets);

        // extract features for them
        objectFeatureExtractor(singlets);
        objectFeatureExtractor(multiplets);

        // save
        singletsSequence.push_back(singlets);
        avgSingletSequence.push_back(AverageObject::average(singlets));
        multipletsSequence.push_back(multiplets);
        
        std::cout << "#T = " << indT 
            << ": #singlets = " << singlets.size()
            << ": #multiplets = " << multiplets.size() << std::endl;

    }

    // generate hypotheses and extract joint features
    hypoSpace(singletsSequence, avgSingletSequence, multipletsSequence);

    return 0;
}
