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

#ifndef __JOINT_FEATURE_EXTRACTOR_HXX__
#define __JOINT_FEATURE_EXTRACTOR_HXX__

#include "Object.hxx"
#include "Event.hxx"
#include "MeasureFactory.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that extracts the joint features using given object features and measures.
 *
 */
class MeasureExtractor {
public:
	/*! Constructor with measures specified 
	 *
	 */
    MeasureExtractor (const std::vector<std::string >& measures)
    {
        measures_ = measures;

        for (int32 ind = 0; ind < measures_.size(); ind ++) {
            MeasureFactory* extractor = MeasureFactory::make(measures_[ind]);
            if (!extractor) {
                std::cerr << "Warning: skipping unknown feature class name " << measures_[ind] << std::endl;
                continue ;
            }

            extractors_.push_back(extractor);
        }
    };
    
	/*! Destructor, delete the joint feature extractors, free the memory
	 *
	 */
    ~MeasureExtractor ()
    {
        for (int32 ind = 0; ind < extractors_.size(); ind ++) 
            delete extractors_[ind];
    };

    /*! Get an extractor by name 
     *  @param measure The extractor (measure) name
     *  @return A pointer to the joint feature extractor
     */
    MeasureFactory* get_extractor(const std::string& measure) 
    {
        std::vector<std::string >::const_iterator it = 
            std::lower_bound(measures_.begin(), measures_.end(), measure);
        
        return get_extractor(static_cast<int32 >(it - measures_.begin()));
    };

    /*! Get an extractor by index 
     *  @param ind The index
     *  @return A pointer to the joint feature extractor
     */
    MeasureFactory* get_extractor(const int32 ind) 
    {
        if (ind < 0)
            return extractors_.front();
        else if (ind >= extractors_.size())
            return extractors_.back();
        
        return extractors_[ind];
    };

    /*! Extract the joint features
     *  @param objects1 The list of objects as the source
     *  @param objects2 The list of objects as the target
     *  @param hypotheses The list of hypothesis
     *  @param e The event definition
     *  @return A Matrix2D object as the joint features
     */
    template<class T1, class T2 >
    Matrix2D operator() (
        const std::vector<T1 >& objects1, 
        const std::vector<T2 >& objects2, 
        const Hypotheses& hypotheses,
        const MeasureSetting& setting) 
    {
        Matrix2D featureMatrix(
            Matrix2D::difference_type(hypotheses.size(), setting.size()), 
            static_cast<MatrixElem >(0));

        for (int32 indH = 0; indH < hypotheses.size(); indH ++) {
            int32 id1 = hypotheses[indH].first;
            int32 id2 = hypotheses[indH].second;
            //std::cout << "extract joint feature for " << id1 << " - " << id2 << std::endl;
            for (int32 indF = 0; indF < setting.size(); indF ++) {
                std::string feature = setting[indF].first;
                std::string measure = setting[indF].second;
                //std::cout << "\tfeature = " << feature << "; measure = " << measure << std::endl;
                MeasureFactory* extractor = get_extractor(measure);
                //std::cout << "\t\textractor" << std::endl;
                featureMatrix(indH, indF) = extractor->operator()(
                    objects1[id1].get_feature(feature), 
                    objects2[id2].get_feature(feature));
            }
        }

        return featureMatrix;
    };

protected:
    std::vector<std::string > measures_;
    std::vector<MeasureFactory* > extractors_;
};
}

#endif /* __JOINT_FEATURE_EXTRACTOR_HXX__ */
