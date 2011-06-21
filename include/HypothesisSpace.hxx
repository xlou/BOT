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

#ifndef __HYPOTHESIS_SPACE_HXX__
#define __HYPOTHESIS_SPACE_HXX__

#include <vector>
#include <map>
#include "FramePair.hxx"
#include "EventConfiguration.hxx"
#include "NearestNeighborGenerator.hxx"
#include "MeasureExtractor.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that creates the entire hypothesis space, including calling 
 *  the hypothesis generation and extracting the joint features.
 *
 */
class HypothesisSpace 
{
public:
	/*! Constructor with the event configuration file specified 
	 *
	 */
    HypothesisSpace(const std::string& confFile) : conf_(confFile)
    {
    };

	/*! Overridden operator(), generate all hypotheses and, for each hypothesis, extract the joint features
     *  @param singlets_vec The sequence of singlets
     *  @param avg_singlet_vec The sequence of average singlet (i.e. the singlet with average features)
     *  @param multiplets_vec The sequence of multiplets
	 */
    void operator()(
        const SingletsSequence& singlets_vec,
        const SingletsSequence& avg_singlet_vec,
        const MultipletsSequence& multiplets_vec)
    {
        //HypothesisGenerator hypoGen(conf_.k(), conf_.d_max());
        MeasureExtractor jfExtractor(conf_.get_measures());
        //ConstraintGenerator constGen;

        // walk through all time steps
        std::cout << "****MeasureExtractor**** " << std::endl;
        for (int32 indT = 0; indT < singlets_vec.size()-1; indT ++) {
            std::cout << "#T = " << indT << " - " << (indT+1) << 
                ": Extract joint features " << std::endl;
            FramePair framepair;
            int32 count1 = singlets_vec[indT].size();
            int32 count2 = singlets_vec[indT+1].size();
            // walk through all events
            for (int32 indE = 0; indE < conf_.count(); indE ++) {
                const std::string& name = conf_.name(indE);
                const EventPairing& pairing = conf_.pairing(indE);
                const MeasureSetting& setting = conf_.setting(indE);
                // generate hypotheses
                Hypotheses hypotheses;
                Matrix2D featureMatrix;
                Matrix2D constraintMatrix;
                if (pairing == std::make_pair<int32, int32 >(1, 1)) {
                    // generate hypotheses
                    hypotheses = generate_hypotheses(singlets_vec[indT], singlets_vec[indT+1]);
                    // generate constraints
                    constraintMatrix = generate_constraint<Singlet, Singlet >(
                        count1, count2, 
                        singlets_vec[indT], singlets_vec[indT+1], hypotheses);
                    // extract joint features
                    featureMatrix = jfExtractor(
                        singlets_vec[indT], singlets_vec[indT+1], 
                        hypotheses, setting);
                }
                else if (pairing == std::make_pair<int32, int32 >(1, 2)) {
                    // generate hypotheses
                    hypotheses = generate_hypotheses(singlets_vec[indT], multiplets_vec[indT+1]);
                    // generate constraints
                    constraintMatrix = generate_constraint<Singlet, Multiplet >(
                        count1, count2, 
                        singlets_vec[indT], multiplets_vec[indT+1], hypotheses);
                    // extract joint features
                    featureMatrix = jfExtractor(
                        singlets_vec[indT], multiplets_vec[indT+1], 
                        hypotheses, setting);
                }
                else if (pairing == std::make_pair<int32, int32 >(2, 1)) {
                    // generate hypotheses
                    hypotheses = generate_hypotheses(multiplets_vec[indT], singlets_vec[indT+1]);
                    // generate constraints
                    constraintMatrix = generate_constraint<Multiplet, Singlet >(
                        count1, count2, 
                        multiplets_vec[indT], singlets_vec[indT+1], hypotheses);
                    // extract joint features
                    featureMatrix = jfExtractor(
                        multiplets_vec[indT], singlets_vec[indT+1], 
                        hypotheses, setting);
                }
                else if (pairing == std::make_pair<int32, int32 >(0, 1)) {
                    // generate hypotheses
                    hypotheses = generate_hypotheses(avg_singlet_vec[indT+1], singlets_vec[indT+1], true);
                    // generate constraints
                    constraintMatrix = generate_constraint<Singlet, Singlet >(
                        count1, count2, 
                        avg_singlet_vec[indT+1], singlets_vec[indT+1], hypotheses);
                    // extract joint features
                    featureMatrix = jfExtractor(
                        avg_singlet_vec[indT+1], singlets_vec[indT+1], 
                        hypotheses, setting);
                }
                else if (pairing == std::make_pair<int32, int32 >(1, 0)) {
                    // generate hypotheses
                    hypotheses = generate_hypotheses(singlets_vec[indT], avg_singlet_vec[indT], true);
                    // generate constraints
                    constraintMatrix = generate_constraint<Singlet, Singlet >(
                        count1, count2, 
                        singlets_vec[indT], avg_singlet_vec[indT], hypotheses);
                    // extract joint features
                    featureMatrix = jfExtractor(
                        singlets_vec[indT], avg_singlet_vec[indT], 
                        hypotheses, setting);
                }
                else {
                    std::cerr << "*Warning* Skipping unsupported event pairing: " << conf_.pairing(indE) << std::endl;
                    continue ;
                }
                framepair.add_event(name, pairing, hypotheses, featureMatrix, constraintMatrix);
            }

            // add this frame pair
            framepairs_.push_back(framepair);
        }

        // preconditioning the features
        precondition();
    };

	/*! Generate the constraint matrix
     *  @param count1 The number of singlets for the source
     *  @param count2 The number of singlets for the target
	 *  @param objects1 A vector of objects as the source
     *  @param objects2 A vector of objects as the target
     *  @param hypotheses A list of hypotheses
     *  @return A Matrix2D object which encodes the constraints
	 */
    template<class T1, class T2 >
    Matrix2D generate_constraint(
        const int32 count1, 
        const int32 count2, 
        const std::vector<T1 >& objects1, 
        const std::vector<T2 >& objects2, 
        const Hypotheses& hypotheses)
    {
        int32 n_variables = hypotheses.size();
        int32 n_constraints = count1 + count2;
        Matrix2D constraintMatrix(
            Matrix2D::difference_type(n_constraints, n_variables), 
            static_cast<MatrixElem >(0));

        for (int32 ind = 0; ind < n_variables; ind ++) {
            const std::vector<int32 >& comp1 = objects1[hypotheses[ind].first].components();
            const std::vector<int32 >& comp2 = objects2[hypotheses[ind].second].components();

            if (comp1[0] != -1) {
                for (int32 ind_id = 0; ind_id < comp1.size(); ind_id ++)
                    constraintMatrix(comp1[ind_id], ind) = 1;
            }
            
            if (comp2[0] != -1) {
                for (int32 ind_id = 0; ind_id < comp2.size(); ind_id ++)
                    constraintMatrix(comp2[ind_id]+count1, ind) = 1;
            }
        }

        return constraintMatrix;
    };


	/*! Generate hypytheses between two lists of objects (singlets or multiplets) 
	 *  @param objects1 A vector of objects as the source
     *  @param objects2 A vector of objects as the target
     *  @param full If true, ignore all neighborhood constraints and create a fully connected bipartite graph
     *  @return A vector of hypothesis, i.e. a pair of indices of objects from those two lists.
	 */
    template<class T1, class T2 >
    Hypotheses generate_hypotheses(
        const std::vector<T1 >& objects1, 
        const std::vector<T2 >& objects2,
        const bool full = false) 
    {
        if (full) {
            Hypotheses hypotheses;
            for (int32 ind1 = 0; ind1 < objects1.size(); ind1 ++)
                for (int32 ind2 = 0; ind2 < objects2.size(); ind2 ++)
                    hypotheses.push_back(std::make_pair(ind1, ind2));

            return hypotheses;
        }

        // search for id pairs in a neighborhood (kNN + distance thresholding)
        NearestNeighborGenerator nn(conf_.k(), conf_.d_max(), false);
        Matrix2D centers1 = VigraSTLInterface::get_centers<T1 >(objects1);
        Matrix2D centers2 = VigraSTLInterface::get_centers<T2 >(objects2);
        return nn(centers1, centers2);
    };

	/*! Return a constant reference to the event configuration 
	 *  @return An EventConfiguration object
	 */
    const EventConfiguration& configuration() const
    {
        return conf_;
    };

	/*! Return a constant reference to the framepairs
	 *  @return An std::vector<FramePair > object
	 */
    const std::vector<FramePair >& framepairs() const
    {
        return framepairs_;
    };

	/*! Return a reference to the framepairs
	 *  @return An std::vector<FramePair > object
	 */
    std::vector<FramePair >& framepairs()
    {
        return framepairs_;
    };

    /*! Precondition the features by normalizing to [0, 1]
     *
	 */
    void precondition()
    {
        // collect all features
        std::vector<Matrix2D > features;
        for (int32 ind = 0; ind < framepairs_.size(); ind ++) {
            const std::vector<Event >& events = framepairs_[ind].events();
            for (int32 indE = 0; indE < events.size(); indE ++) {
                if (ind == 0)
                    features.push_back(events[indE].features());
                else
                    features[indE] = joinVertically(features[indE], events[indE].features());
            }
        }

        // compute their 5% and 95% quantile
        std::vector<Matrix2D > vec_lb;
        std::vector<Matrix2D > vec_ub;
        for (int32 ind = 0; ind < features.size(); ind ++) {
            int32 n_rows = features[ind].shape(0);
            int32 n_cols = features[ind].shape(1);
            int32 nth = static_cast<int32 >(n_rows * 0.10);

            Matrix2D lb(Shape2D(1, n_cols));
            Matrix2D ub(Shape2D(1, n_cols));
            for (int32 col = 0; col < n_cols; col ++) {
                Shape2D p(0, col);
                Shape2D q(n_rows, col+1);
                Matrix2D tmp(features[ind].subarray(p, q));

                std::nth_element(tmp.begin(), tmp.begin() + nth, tmp.end());
                lb[col] = tmp[nth];

                std::nth_element(tmp.begin(), tmp.begin() + nth, tmp.end(), std::greater<MatrixElem >());
                ub[col] = tmp[nth];
            }
            vec_lb.push_back(lb);
            vec_ub.push_back(ub);
        }

        // normalize
        for (int32 ind = 0; ind < framepairs_.size(); ind ++) {
            std::vector<Event >& events = framepairs_[ind].events();
            for (int32 indE = 0; indE < events.size(); indE ++) {
                Matrix2D& F = events[indE].features();
                Matrix2D mat_lb = 
                    repeatMatrix(vec_lb[indE], F.shape(0), 1);
                Matrix2D mat_ub = 
                    repeatMatrix(vec_ub[indE], F.shape(0), 1) +
                    1e-8;   // for numerical stability reasons
                Matrix2D X = pdiv(static_cast<MatrixElem >(2) * F - mat_lb - mat_ub, mat_ub - mat_lb);
                F = static_cast<MatrixElem >(1) / 
                    (static_cast<MatrixElem >(1) + exp(-X)) - static_cast<MatrixElem >(1);
            }
        }

    };
protected:
    std::vector<FramePair > framepairs_;
    EventConfiguration conf_;
};

}

#endif /* __HYPOTHESIS_SPACE_HXX__ */
