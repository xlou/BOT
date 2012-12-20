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
#include "vigra/hdf5impex.hxx"
#include "SolutionCoder.hxx"
#include "TrainingData.hxx"
#include "TrackingTrainer.hxx"
#include "HDF5ReaderWriter.hxx"

using namespace bot;

int main(int argc, char* argv[])
{
        if (argc != 3) {
                std::cerr << "********Input error, please use********" << std::endl <<
                "                  ./TrackPredictor [hdf5 data file] [ini configuration file]" << std::endl;
		return EXIT_FAILURE;
        }

    std::string filename(argv[1]);
    // load the image sequence
    std::vector<Matrix2D > images, segmentations;
    TrainingData training;
    HDF5ReaderWriter::load(filename, images, segmentations);
    std::cout << "****Loading the images/segmentations****" << std::endl;
    HDF5ReaderWriter::load(filename, training);
    std::cout << "****Loading the training data****" << std::endl;

    // get the context
    Context context(images);
    std::cout << "****Computing the Context****" << std::endl << context << std::endl << std::endl;

    // load the configuration
    HypothesisSpace space(argv[2]);
    EventConfiguration conf = space.configuration();

    // create singlets/muliplets and extract object features
    std::cout << "****Extracting singlets and multiplets****" << std::endl;
    SingletsSequence singlets_vec;
    SingletsSequence avg_singlet_vec;
    MultipletsSequence multiplets_vec;
    SingletsGenerator singletsGenerator;
    MultipletsGenerator multipletsGenerator(conf.k(), conf.d_max());
    ObjectFeatureExtractor extractor(conf.get_feature_names(), context);
    for (int32 indT = 0; indT < images.size(); indT ++) {
        // generate singlets and multiplets
        Singlets singlets = singletsGenerator(images[indT], segmentations[indT]);
        Multiplets multiplets = multipletsGenerator(images[indT], segmentations[indT], singlets);

        // extract features for them
        extractor(singlets);
        extractor(multiplets);
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
    const std::vector<FramePair >& framepairs = space.framepairs();

    // parse the training data
    std::cout << "****Parsing the training data****" << std::endl;
    SolutionCoder coder;
    int32 nTr = training.times().size();
    for (int32 ind = 0; ind < nTr; ind ++) {
        int32 time = training.times()[ind];
        std::cout << "****time = " << time << "****" << std::endl;
        const LabelAssociations& association = training.associations()[ind];

        const std::vector<Event >& events = framepairs[time].events();
        const Singlets& singlets1 = singlets_vec[time];
        const Singlets& singlets2 = singlets_vec[time+1];
        const Multiplets& multiplets1 = multiplets_vec[time];
        const Multiplets& multiplets2 = multiplets_vec[time+1];

        Solution solution;
        coder.decode(
            association, 
            events, 
            singlets1, singlets2,
            multiplets1, multiplets2,
            solution);
        training.solutions().push_back(solution);
    }

    // start the training
    TrackingTrainer trainer;
    const std::vector<Matrix2D > null_vector;
    std::vector<Matrix2D > weights = conf.weights(0.5);
    std::string msg = trainer(training, framepairs, weights, true);
    std::cout << "Training returns: " << msg << std::endl;
    conf.weights() = weights;

    // printe intermediate results: weights, epsilons, losses
    std::cout << "************Intermediate output: **************" << std::endl;
    std::cout << "Weights: " << std::endl << trainer.weights() << std::endl << std::endl;
    std::cout << "Epsilons: " << std::endl << trainer.epsilons() << std::endl << std::endl;
    std::cout << "Losses: " << std::endl << trainer.losses() << std::endl << std::endl;

    // print the final weights
    std::cout << "Learned weights: " << std::endl;
    conf.print();


    return EXIT_SUCCESS;
}
