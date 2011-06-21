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

#ifndef __OBJECT_FEATURE_CENTERS_HXX__
#define __OBJECT_FEATURE_CENTERS_HXX__

#include "objectFeatures.hxx"
#include "ObjectFeatureFactory.hxx"
#include "InputOutput.hxx"

namespace bot {

/*! A class for extracting the centers for each label
 *  This class extends class ObjectFeatureFactory and overrides its virtual function extract() and shape()
 */
class ObjectFeatureCenters : public ObjectFeatureFactory {
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() {
        return std::string("ObjectFeatureCenters");
    }

    /*! A static function that returns the centers for each label
     *  @return A Matrix2D::difference_type variable that represents the shape of the feature matrix
     */
    Matrix2D::difference_type shape(int dim = 2) {
        return Matrix2D::difference_type(2, dim);
    }

	/*! Implementation of the virtual function
     *  @param feature_mat The object feature to be returned
     *  @param obj The input object
     *  @param context The global context
	 */
    void extract(Matrix2D& feature_mat, const Object& obj, const Context& context = Context()) {
        if (obj.is_singlet()) {
            feature_mat.reshape(Matrix2D::difference_type(1, obj.dim()), 0);
            rowVector(feature_mat, 0) = obj.center();
            return ;
        }

        feature_mat.reshape(shape(obj.dim()), 0);
        std::map<MatrixElem, Matrix2D > label_to_centers;
        std::map<MatrixElem, MatrixElem > label_to_counts;
        std::vector<MatrixElem > labels;
        for (int32 ind = 0; ind < obj.labels().size(); ind ++) {
            MatrixElem label = obj.labels()[ind];
            if (label == 0)
                continue ;

            if (label_to_centers.find(label) == label_to_centers.end()) {
                label_to_centers[label] = Matrix2D(
                    Matrix2D::difference_type(1, obj.dim()), 
                    static_cast<MatrixElem >(0));
                label_to_counts[label] = 0;
                labels.push_back(label);
            }

            label_to_centers[label] = label_to_centers[label] + rowVector(obj.pixels(), ind);
            label_to_counts[label] += 1;
        }

        rowVector(feature_mat, 0) = label_to_centers[labels[0]] / 
            static_cast<MatrixElem >(label_to_counts[labels[0]]);
        rowVector(feature_mat, 1) = label_to_centers[labels[1]] / 
            static_cast<MatrixElem >(label_to_counts[labels[1]]);
    }
};

}
#endif /* __OBJECT_FEATURE_CENTERS_HXX__ */
