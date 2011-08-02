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
#include "CPLEXSolverSystem.hxx"
#include "TrackingPredictor.hxx"
#include "HDF5ReaderWriter.hxx"

using namespace bot;

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
    HypothesisSpace space("../data/event-configuration-cell.ini");
    EventConfiguration conf = space.configuration();

    // create singlets/muliplets and extract object features
    std::cout << "****Extracting singlets and multiplets" << std::endl;
    SingletsSequence singlets_vec;
    SingletsSequence avg_singlet_vec;
    MultipletsSequence multiplets_vec;
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
        singlets_vec.push_back(singlets);
        avg_singlet_vec.push_back(AverageObject::average(singlets));
        multiplets_vec.push_back(multiplets);
        
        std::cout << "#T = " << indT 
            << ": #singlets = " << singlets.size()
            << ": #multiplets = " << multiplets.size() << std::endl;

    }

    // generate hypotheses and extract joint features
    space(singlets_vec, avg_singlet_vec, multiplets_vec);


    /*
     * perform tracking
     */
    // create cplex solver
    CPLEXSolverSystem solver;
    TrackingPredictor predictor;
    std::vector<FramePair >& framepairs = space.framepairs();
    const std::vector<Matrix2D > null_vector;
    for (int32 ind = 0; ind < framepairs.size(); ind ++) {
        Solution& solution = framepairs[ind].solution(); 
        std::string msg = predictor(framepairs[ind], conf.weights(), solution, null_vector);
        std::cout << "T = " << ind << "; cplex returns: " << msg << std::endl;
    }

    // save the solutions to the hdf file
    HDF5ReaderWriter::save(filename, framepairs, singlets_vec, multiplets_vec);

    return 0;
}
