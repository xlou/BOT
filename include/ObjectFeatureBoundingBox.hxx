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

#ifndef __OBJECT_FEATURE_BOUNDING_BOX_HXX__
#define __OBJECT_FEATURE_BOUNDING_BOX_HXX__

#include "objectFeatures.hxx"
#include "ObjectFeatureFactory.hxx"

namespace bot {

/*! A class for extracting the bounding box
 *  This class extends class ObjectFeatureFactory and overrides its virtual function extract() and shape()
 */
class ObjectFeatureBoundingBox : public ObjectFeatureFactory {
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() {
        return std::string("ObjectFeatureBoundingBox");
    }

    /*! A static function that returns the shape of the feature matrix
     *  @return A Matrix2D::difference_type variable that represents the shape of the feature matrix
     */
    Matrix2D::difference_type shape(int dim = 2) {
        return Matrix2D::difference_type(1, 2*dim);
    }

	/*! Implementation of the virtual function
     *  @param feature_mat The object feature to be returned
     *  @param obj The input object
     *  @param context The global context
	 */
    void extract(Matrix2D& feature_mat, const Object& obj, const Context& context = Context()) {
        feature_mat.reshape(shape(obj.dim()), 0);

        three_set coords = VigraSTLInterface::pixels2three_set(obj.pixels());
        value_set intens = VigraSTLInterface::values2value_set(obj.values());
        features::ObjectBoundingBox<unsigned short, 3 > o(coords, intens);
        features::features_type ret = o.get();

        if (obj.dim() == 3)
            std::copy (ret.begin(), ret.begin()+2*obj.dim(), feature_mat.begin());
        else {
            feature_mat(0, 0) = ret[0];
            feature_mat(0, 1) = ret[1];
            feature_mat(0, 2) = ret[3];
            feature_mat(0, 3) = ret[4];
        }
    }
};

}
#endif /* __OBJECT_FEATURE_BOUNDING_BOX_HXX__ */
