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

#ifndef __JOINT_FEATURE_FACTORY_HXX__
#define __JOINT_FEATURE_FACTORY_HXX__

#include "TypeDefinition.hxx"

namespace bot {

/*! A class as the base class for joint feature extractor
 *
 */
class MeasureFactory
{
public:
	/*! Default constructor 
	 *
	 */
    MeasureFactory()
    {
        // call the default initialization, e.g. setting default parameters
        initialize();
    };

	/*! A virtual function that extract the feature
     *  @param feature1 The object feature of the source
     *  @param feature2 The object feature of the target
	 *  @return A MatrixElem value as the joint feature
	 */
    virtual MatrixElem operator()(const Matrix2D& feature1, const Matrix2D& feature2) = 0;

    /*! A virtual function that initialize the class
     *
	 */
    virtual void initialize() {};
    
	/*! A static function that implements the factory method with a given class name
     *  @param classname The name of the joint feature extractor
	 *  @return A point to a newly created joint feature extractor
	 */
    static MeasureFactory *make(const std::string& classname);

    /*! A static function that implements the factory method with a given class name
     *  @param classname The name of the joint feature extractor
     *  @param param The parameters
	 *  @return A point to a newly created joint feature extractor
	 */
    static MeasureFactory *make(const std::string& classname, const Matrix2D& param);

    /*! Return a constant reference to the parameters
	 *  @return A Matrix2D object as the parameter matrix
	 */
    const Matrix2D& param() const 
    {
        return param_;
    };

     /*! Set the parameters
	 *  @param param The parameter matrix to be set
	 */
    void set_param(const Matrix2D& param) 
    {
        param_ = param;
    };
protected:
    Matrix2D param_;
};

}

#endif /* __JOINT_FEATURE_FACTORY_HXX__ */
