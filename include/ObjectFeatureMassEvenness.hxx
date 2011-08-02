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

#ifndef __OBJECT_FEATURE_MASS_EVENNESS_HXX__
#define __OBJECT_FEATURE_MASS_EVENNESS_HXX__

#include "vigra/matrix.hxx"
#include "ObjectFeatureFactory.hxx"

using namespace vigra::linalg;

namespace bot {

/*! A class for extracting the mass eveness for a multiplet
 *  This class extends class ObjectFeatureFactory and overrides its virtual function extract() and shape()
 */
class ObjectFeatureMassEvenness : public ObjectFeatureFactory {
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() {
        return std::string("ObjectFeatureMassEvenness");
    }

    /*! A static function that returns the shape of the return feature
     *  @return A Matrix2D::difference_type variable that represents the shape of the return feature
     */
    Matrix2D::difference_type shape(int dim = 2) {
        return Matrix2D::difference_type(1, 1);
    }

	/*! Implementation of the virtual function
     *  @param feature_mat The object feature to be returned
     *  @param obj The input object
     *  @param context The global context
	 */
    void extract(Matrix2D& feature_mat, const Object& obj, const Context& context) {
        feature_mat.reshape(shape(obj.dim()), 0);

        MatrixElem mass_evenness;
        int num_zeros = std::count_if(
            obj.labels().begin(), 
            obj.labels().end(), 
            std::bind2nd(std::equal_to<MatrixElem >(), static_cast<MatrixElem >(0)));
        if (num_zeros == 0) {
            mass_evenness = static_cast<MatrixElem >(0.0);
        }
        else {
            Matrix2D NZ(obj.labels());   // indicator matrix of non-zeros
            std::replace_if(
                NZ.begin(), 
                NZ.end(), 
                std::bind2nd(std::not_equal_to<MatrixElem >(), static_cast<MatrixElem >(0)), 
                static_cast<MatrixElem >(1));

            Matrix2D Z = static_cast<MatrixElem >(1) - NZ;
            MatrixElem mean_inten_z = dot(Z, obj.values())/num_zeros;
            MatrixElem mean_inten_nz = dot(NZ, obj.values())/(obj.count() - num_zeros);
            mass_evenness = static_cast<MatrixElem >(mean_inten_z) / 
                static_cast<MatrixElem >(mean_inten_nz);
        }
        
        feature_mat(0, 0) = mass_evenness;
    }
};

}
#endif /* __OBJECT_FEATURE_MASS_EVENNESS_HXX__ */
