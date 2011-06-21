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

#ifndef __JOINT_FEATURE_APPOINTMENT_HXX__
#define __JOINT_FEATURE_APPOINTMENT_HXX__

#include <functional>
#include <numeric>
#include "vigra/matrix.hxx"
#include "MeasureFactory.hxx"

using namespace vigra::linalg;

namespace bot {

/*! A class for directly appointing the first feature to the output feature
 *  This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
 */
class MeasureAppointmentLeft : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureAppointmentLeft");
    };

    /*! Directly assign the value of the first (left) feature as the output
     *  @param feature1 a bot::Matrix2D type feature vector (matrix)
     *  @param feature2 a bot::Matrix2D type feature vector (matrix)
     *  @return A MatrixElem value
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        return feature1[0];
    };
};

/*! A class for directly appointing the second feature to the output feature
 *  This class extends class MeasureFactory and overrides its virtual function operator() and initialize()
 */
class MeasureAppointmentRight : public MeasureFactory 
{
public:
    /*! A static function that returns the class name
     *  @return A std::string variable that represents the class name
     */
    static std::string getClassName() 
    {
        return std::string("MeasureAppointmentRight");
    };

    /*! Directly assign the value of the second (right) feature as the output
     *  @param feature1 a bot::Matrix2D type feature vector (matrix)
     *  @param feature2 a bot::Matrix2D type feature vector (matrix)
     *  @return A MatrixElem value
     */
    MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) 
    {
        return feature2[0];
    };
};

}
#endif /* __JOINT_FEATURE_APPOINTMENT_HXX__ */
