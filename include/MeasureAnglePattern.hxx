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

#ifndef __JOINT_FEATURE_ANGLE_PATTERN_HXX__
#define __JOINT_FEATURE_ANGLE_PATTERN_HXX__

#include <functional>
#include <numeric>
#include "vigra/matrix.hxx"
#include "MeasureFactory.hxx"

using namespace vigra::linalg;

namespace bot {

/*! A class for computing the angle between one object and the other two objects.
 *  This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
 */
class MeasureAnglePattern : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureAnglePattern");
    };

    /*! Initialize this class, e.g. setting the default parameters
     *
     */
    void initialize()
    {
    };

    /*! Return the the angle between one object and the other two objects
     *  @return A MatrixElem variable
     */
    MatrixElem dot_product(const Matrix2D& center, const Matrix2D& centers)
    {
        Matrix2D v1 = center - rowVector(centers, 0);
        MatrixElem norm1 = norm(v1);
        Matrix2D v2 = center - rowVector(centers, 1);
        MatrixElem norm2 = norm(v2);

        if (norm1 == 0 || norm2 == 0) 
            return 0;
        
        return dot(v1 /= norm1, v2 /= norm2);
    };

    /*! Compute the doc product of the vectors between one object center and the other two object centers
     *  @param feature1 a bot::Matrix2D type feature vector (matrix)
     *  @param feature2 a bot::Matrix2D type feature vector (matrix)
     *  @return A MatrixElem value as the dot product
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        // must be between one center and the other two centers
        if (!(feature1.shape(0) == 2 && feature2.shape(0) == 1) &&
            !(feature1.shape(0) == 1 && feature2.shape(0) == 2))
            return 0;

        MatrixElem d;
        if (feature1.shape(0) == 1)
            d = dot_product(feature1, feature2);
        else
            d = dot_product(feature2, feature1);

        return d;
    };
};

}
#endif /* __JOINT_FEATURE_ANGLE_PATTERN_HXX__ */
