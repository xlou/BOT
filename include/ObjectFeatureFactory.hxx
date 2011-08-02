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

#ifndef __OBJECT_FEATURE_FACTORY_HXX__
#define __OBJECT_FEATURE_FACTORY_HXX__

#include <iostream>
#include "Object.hxx"
#include "Context.hxx"
#include "TypeDefinition.hxx"

namespace bot {

/*! A class as the base class for object feature extractor
 *
 */
class ObjectFeatureFactory
{
public:
	/*! Default constructor 
	 *
	 */
    ObjectFeatureFactory() {};

	/*! A virtual function that extract the object feature
     *  @param feature_mat The object feature to be returned
     *  @param obj The input object
     *  @param context The global context
	 */
    virtual void extract(Matrix2D& feature_mat, const Object& obj, const Context& context) = 0;

	/*! Return the shape (size) of the feature matrix
     *  @param dim The input data dimension
     *  @return A Matrix2D::difference_type object as the shape
	 */
    virtual Matrix2D::difference_type shape(int dim = 2) = 0;
    
	/*! A static function that implements the factory method with a given class name
     *  @param classname The name of the object feature extractor
	 *  @return A point to a newly created object feature extractor
	 */
    static ObjectFeatureFactory *make(const std::string& classname);

    /*! A static function that implements the factory method with a given class name
     *  @param classname The name of the object feature extractor
     *  @param param The parameters
	 *  @return A point to a newly created object feature extractor
	 */
    static ObjectFeatureFactory *make(const std::string& classname, const Matrix2D& param);

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

#endif /* __OBJECT_FEATURE_FACTORY_HXX__ */
