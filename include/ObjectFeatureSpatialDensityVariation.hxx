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

#ifndef __OBJECT_FEATURE_SPATIAL_DENSITY_VARIATION_HXX__
#define __OBJECT_FEATURE_SPATIAL_DENSITY_VARIATION_HXX__

#include <vigra/matrix.hxx>
#include <algorithm>
#include <functional>
#include "ObjectFeatureFactory.hxx"

using namespace vigra::linalg;

namespace bot {

/*! A class for extracting the spatial density variation
 *  This class extends class ObjectFeatureFactory and overrides its virtual function extract() and shape()
 */
class ObjectFeatureSpatialDensityVariation : public ObjectFeatureFactory {
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() {
        return std::string("ObjectFeatureSpatialDensityVariation");
    }

    /*! A static function that returns the shape of the feature matrix
     *  @return A Matrix2D::difference_type variable that represents the shape of the feature matrix
     */
    Matrix2D::difference_type shape(int dim = 2) {
        return Matrix2D::difference_type(1, 1);
    }

	/*! Implementation of the virtual function
     *  @param feature_mat The object feature to be returned
     *  @param obj The input object
     *  @param context The global context
	 */
    void extract(Matrix2D& feature_mat, const Object& obj, const Context& context = Context()) {
        const Matrix2D& pixels = obj.pixels();
        const Matrix2D& values = obj.values();

        // convert intensity to probability
        Matrix2D p(values.shape());
        sdiv(values, std::accumulate(values.begin(), values.end(), 0.0), p);

        // get mass center
        Matrix2D X = pmul(pixels, repeatMatrix(p, 1, obj.dim()));
        Matrix2D mass_center(Matrix2D::difference_type(1, obj.dim()));
        columnStatistics(X *= obj.count(), mass_center);
        X = pixels - repeatMatrix(mass_center, obj.count(), 1);

        // computer points' distance to the mass center
        Matrix2D row_mean(values.shape()), row_std(values.shape()), row_norm(values.shape());
        rowStatistics(X, row_mean, row_std, row_norm);
        X = pmul(row_norm, p);

        // return the feature value
        feature_mat.reshape(shape(obj.dim()), 0);
        feature_mat[0] = MatrixElem(std::accumulate(X.begin(), X.end(), 0.0)) / 
            MatrixElem(std::accumulate(values.begin(), values.end(), 0.0));
    }
};

}
#endif /* __OBJECT_FEATURE_SPATIAL_DENSITY_VARIATION_HXX__ */
