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

#ifndef __JOINT_FEATURE_L1L2_DISTANCE_HXX__
#define __JOINT_FEATURE_L1L2_DISTANCE_HXX__

#include <numeric>
#include "vigra/matrix.hxx"
#include "MeasureFactory.hxx"

using namespace vigra::linalg;

namespace bot {

/*! Return the norm (l1 or l2) of a given matrix 
 *  @param feature The input matrix
 *  @param p The p-norm
 *  @return An MatrixElem value as the norm
 */
MatrixElem get_norm(const Matrix2D& feature, const MatrixElem p = 2)
{
    MatrixElem d;
    if (p == 1) {
        Matrix2D X = abs(feature);
        d = std::accumulate(X.begin(), X.end(), static_cast<MatrixElem >(0));
    }
    else {
        d = norm(feature);
    }
    
    return d;
};

/*! Return the distance (l1 or l2) of two features
 *  @param feature1 The feature as the source
 *  @param feature2 The feature as the target
 *  @param p The p-norm
 *  @param normalize Whether the result shall be normalized
 *  @return An MatrixElem value as the distance
 */
MatrixElem distance(const Matrix2D& feature1, const Matrix2D& feature2, 
    const MatrixElem p, const MatrixElem normalize) 
{
    MatrixElem d = get_norm(feature1 - feature2, p);
    if (normalize != 0) {
        MatrixElem norm1 = get_norm(feature1, p);
        if (norm1 == 0 && d == 0)  
            return 0;   // the two features much be both zero

        MatrixElem norm2 = get_norm(feature2, p);
        d = d / std::max(norm1, norm2);
    }

    return d;
};

/*! A class for computing the Euclidean distance
    This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
*/
class MeasureEuclideanDistance : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureEuclideanDistance");
    };

    /*! Compute the Euclidean distance of given two feature vectors (or matrics)
     *  @param feature1 The feature as the source
     *  @param feature2 The feature as the target
     *  @return An MatrixElem value as the distance
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        return distance(feature1, feature2, 2, 0);
    };
};

/*! A class for computing the normalized Euclidean distance
    This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
*/
class MeasureNormalizedEuclideanDistance : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureNormalizedEuclideanDistance");
    };

    /*! Compute the normalized Euclidean distance of given two feature vectors (or matrics)
     *  @param feature1 The feature as the source
     *  @param feature2 The feature as the target
     *  @return An MatrixElem value as the distance
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        return distance(feature1, feature2, 2, 1);
    };
};

/*! A class for computing the Manhattan distance
    This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
*/
class MeasureManhattanDistance : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureManhattanDistance");
    };

    /*! Compute the Manhattan distance of given two feature vectors (or matrics)
     *  @param feature1 The feature as the source
     *  @param feature2 The feature as the target
     *  @return An MatrixElem value as the distance
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        return distance(feature1, feature2, 1, 0);
    };
};

/*! A class for computing the normalized Manhattan distance
    This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
*/
class MeasureNormalizedManhattanDistance : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureNormalizedManhattanDistance");
    };

    /*! Compute the normalized Manhattan distance of given two feature vectors (or matrics)
     *  @param feature1 The feature as the source
     *  @param feature2 The feature as the target
     *  @return An MatrixElem value as the distance
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        return distance(feature1, feature2, 1, 1);
    };
};

}
#endif /* __JOINT_FEATURE_L1L2_DISTANCE_HXX__ */
