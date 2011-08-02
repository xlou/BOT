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

#ifndef __TRAINING_DATA_HXX__
#define __TRAINING_DATA_HXX__

#include <vector>
#include <vigra/matrix.hxx>
#include "TypeDefinition.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! This class manages the trainig data
 *
 */
class TrainingData 
{
public:
    /*! Add a new training sample
     *  @param time The index of the left frame
     *  @param associations A LabelAssociations object 
     */
    void add(const int32 time, const LabelAssociations& associations)
    {
        times_.push_back(time);
        associations_.push_back(associations);
    };

    /*! Return a constant reference to the times
     *  @return A constant reference
     */
    const std::vector<int32 >& times() const
    {
        return times_;
    };

    /*! Return a constant reference to the solutions
     *  @return A constant reference
     */
    const std::vector<Solution >& solutions() const
    {
        return solutions_;
    };
    
    /*! Return a reference to the solutions
     *  @return A reference
     */
    std::vector<Solution >& solutions()
    {
        return solutions_;
    };

    /*! Return a constant reference to the associations
     *  @return A constant reference
     */
    const std::vector<LabelAssociations >& associations() const
    {
        return associations_;
    };

protected:
    std::vector<int32 > times_;
    std::vector<Solution > solutions_;
    std::vector<LabelAssociations > associations_;
};

}

#endif /* __TRAINING_DATA_HXX__ */
