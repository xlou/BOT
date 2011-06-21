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

#ifndef __OBJECT_FEATURE_ECCENTRICITY_HXX__
#define __OBJECT_FEATURE_ECCENTRICITY_HXX__

#include "objectFeatures.hxx"
#include "ObjectFeatureFactory.hxx"

namespace bot {

/*! \brief ObjectFeatureEccentricity: A class for extracting the eccentricity of the shape
    This class extends class ObjectFeatureFactory and overrides its virtual function extract() and shape()
*/
class ObjectFeatureEccentricity : public ObjectFeatureFactory {
public:
    /*! A static function that returns the class name
        \return A std::string variable that represents the class name
    */
    static std::string getClassName() {
        return std::string("ObjectFeatureEccentricity");
    };

    /*! A static function that returns the shape of the feature matrix
        \return A FeatureMatrix::difference_type variable that represents the shape of the feature matrix
    */
    FeatureMatrix::difference_type shape(int dim = 2) {
        return FeatureMatrix::difference_type(1, 1);
    };

    /*! Extract the eccentricity of a given object
        \param obj a bot::Object type type object
        \param obj a bot::FeatureMatrix type eccentricity of the shape (return value)
    */
    void extract(FeatureMatrix& feature_mat, const Object& obj, const Context& context = Context()) {
        feature_mat.reshape(shape(obj.dim()), 0);

        three_set coords = VigraSTLInterface::pixels2three_set(obj.pixels());
        value_set intens = VigraSTLInterface::values2value_set(obj.values());
        features::ObjectPrincipalComponents<unsigned short, 3 > o(coords, intens);
        features::features_type ret = o.get();

        feature_mat(0, 0) = ret[0]/ret[obj.dim()-1];
    };
};

}
#endif /* __OBJECT_FEATURE_ECCENTRICITY_HXX__ */
