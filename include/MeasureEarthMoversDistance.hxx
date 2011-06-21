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

#ifndef __JOINT_FEATURE_EARTH_MOVERS_DISTANCE_HXX__
#define __JOINT_FEATURE_EARTH_MOVERS_DISTANCE_HXX__

#include <functional>
#include <numeric>
#include "MeasureFactory.hxx"
#include "VigraSTLInterface.hxx"
#include "FastEMD/emd_hat.hpp"

using namespace vigra::linalg;

namespace bot {

/*! A class for computing the earth mover's distance
 *  This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
 */
class MeasureEarthMoversDistance : public MeasureFactory 
{
public:
	/*! Default constructor, initialze default parameters 
	 *
	 */
    MeasureEarthMoversDistance()
    {
        // set the default parameters
        Matrix2D param(Matrix2D::difference_type(1, 1), static_cast<MatrixElem >(0));
        param[0] = 1e10;       // distance threshold is 10
        set_param(param);
    };

    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureEarthMoversDistance");
    };

    /*! Return the distance threshold
     *  @return The distance threshold
     */
    MatrixElem threshold()
    {
        return param_[0];
    }

    /*! Update the thresholded ground distance matrix
     *  @param n_bins The number of bins
     */
    void update_distance_matrix(int n_bins)
    {
        C_.clear();
        for (MatrixElem ind1 = 0; ind1 < n_bins; ind1 ++) {
            std::vector<MatrixElem > row;
            for (MatrixElem ind2 = 0; ind2 < n_bins; ind2 ++) 
                //row.push_back(std::min(threshold(), std::abs(ind2 - ind1)));
                row.push_back(std::abs(ind2 - ind1));

            C_.push_back(row);
        }
    };

    /*! Compute the earth mover's distance of given two feature vectors (or matrics)
     *  @param feature1 a bot::Matrix2D type feature vector (matrix)
     *  @param feature2 a bot::Matrix2D type feature vector (matrix)
     *  @return A MatrixElem value
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        std::vector<MatrixElem > hist1 = VigraSTLInterface::view_to_vector(feature1);
        std::vector<MatrixElem > hist2 = VigraSTLInterface::view_to_vector(feature2);

        if (C_.size() != hist1.size()) 
            update_distance_matrix(hist1.size());

        MatrixElem emd = emd_hat_gd_metric<MatrixElem >()(hist1, hist2, C_);

        return emd;
    };

private:
    std::vector<std::vector<MatrixElem > > C_;
};

}
#endif /* __JOINT_FEATURE_EARTH_MOVERS_DISTANCE_HXX__ */
