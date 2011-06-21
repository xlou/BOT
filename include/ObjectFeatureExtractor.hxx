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

#ifndef __OBJECT_FEATURE_EXTRACTOR_HXX__
#define __OBJECT_FEATURE_EXTRACTOR_HXX__

#include "Object.hxx"
#include "Event.hxx"
#include "Context.hxx"
#include "ObjectFeatureFactory.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that extracts all required object features for object
 *
 */
class ObjectFeatureExtractor {
public:
	/*! Constructor with object feature names and the global context specified 
	 *  @param names The object feature names
     *  @parma context The global context
	 */
    ObjectFeatureExtractor (const std::vector<std::string >& names, const Context& context)
    {
        names_ = names;

        for (int32 ind = 0; ind < names_.size(); ind ++) {
            ObjectFeatureFactory* extractor = ObjectFeatureFactory::make(names_[ind]);
            if (!extractor) {
                std::cerr << "*Warning* Skipping unknown feature class name " << names_[ind] << std::endl;
                continue ;
            }

            extractors_.push_back(extractor);
        }

        context_ = context;
    };
    
	/*! Destructor, delete the feature extractors and free the memory
	 *
	 */
    ~ObjectFeatureExtractor ()
    {
        for (int32 ind = 0; ind < extractors_.size(); ind ++) 
            delete extractors_[ind];
    };

	/*! Overridden operator(), extract features for an object
     *  @param object The input object
	 */
    template<class T >
    void operator() (T& object) 
    {
        for (int32 ind = 0; ind < extractors_.size(); ind ++) {
            Matrix2D featureMatrix;
            extractors_[ind]->extract(featureMatrix, object, context());
            object.add_feature(names()[ind], featureMatrix);
        }
    };

    /*! Overridden operator(), extract features for a vector of objects
     *  @param object The vector of objects
	 */
    template<class T >
    void operator() (std::vector<T >& objects) 
    {
        for (int32 ind = 0; ind < objects.size(); ind ++) 
            operator()(objects[ind]);
    };

	/*! Return a constant reference to the object feature names
	 *  @return A vector of std::string as the object feature names
	 */
    const std::vector<std::string >& names() const
    {
        return names_;
    };

	/*! Return a constant reference to the global context
	 *  @return A Context object as the global context
	 */
    const Context& context() const
    {
        return context_;
    };
    
protected:
    Context context_;
    std::vector<std::string > names_;
    std::vector<ObjectFeatureFactory* > extractors_;
};

}

#endif /* __OBJECT_FEATURE_EXTRACTOR_HXX__ */
