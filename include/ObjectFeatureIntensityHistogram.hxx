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

#ifndef __OBJECT_FEATURE_INTENSITY_HISTOGRAM_HXX__
#define __OBJECT_FEATURE_INTENSITY_HISTOGRAM_HXX__

#include <math.h>
#include <algorithm>
#include "objectFeatures.hxx"
#include "ObjectFeatureFactory.hxx"

namespace bot {

/*! A class for extracting the intensity histogram
 *  This class extends class ObjectFeatureFactory and overrides its virtual function extract() and shape()
 */
class ObjectFeatureIntensityHistogram : public ObjectFeatureFactory {
public:
    /*! Overridden constructor: initialize the parameters
     *  @return A bot::ObjectFeatureIntensityHistogram class instance
     */
    ObjectFeatureIntensityHistogram() {
        param_.reshape(Matrix2D::difference_type(1, 3));
        param_[0] = 0;
        param_[1] = 5;
        param_[2] = 255;
    }

    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() {
        return std::string("ObjectFeatureIntensityHistogram");
    }

    /*! A static function that returns the shape of the feature matrix
     *  @return A Matrix2D::difference_type variable that represents the shape of the feature matrix
     */
    Matrix2D::difference_type shape(int dim = 2) {
        int no_bins = static_cast<int >(1 + ceil(param_[2] - param_[0])/param_[1]);
        return Matrix2D::difference_type(1, static_cast<MatrixElem >(no_bins));
    }

	/*! Implementation of the virtual function
     *  @param feature_mat The object feature to be returned
     *  @param obj The input object
     *  @param context The global context
	 */
    void extract(Matrix2D& feature_mat, const Object& obj, const Context& context = Context()) {
        int no_bins = 101;
        MatrixElem interval = ceil((context.max() - context.min()) / static_cast<MatrixElem >(no_bins)) + 1;
        //feature_mat.reshape(shape(obj.dim()), 0);
        feature_mat.reshape(
            Matrix2D::difference_type(1, no_bins), 
            static_cast<MatrixElem >(0));

        for (int32 ind = 0; ind < obj.count(); ind ++) {
            MatrixElem v = static_cast<MatrixElem >(obj.values()[ind] - context.min());
            int32 bin = std::min(no_bins - 1, static_cast<int32 >(v/interval));
            feature_mat[bin] = feature_mat[bin] + 1;
        }

        MatrixElem sum_ = std::accumulate(
            feature_mat.begin(), 
            feature_mat.end(),
            static_cast<MatrixElem >(0));
        feature_mat /= sum_;
    }
};

}
#endif /* __OBJECT_FEATURE_INTENSITY_HISTOGRAM_HXX__ */
