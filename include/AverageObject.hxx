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

#ifndef __AVERAGE_OBJECT_HXX__
#define __AVERAGE_OBJECT_HXX__

#include <vector>
#include <vigra/matrix.hxx>
#include "Object.hxx"
#include "TypeDefinition.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that computes the average features for a given list of objects.
 *
 */
class AverageObject 
{
public:
    /*! Compute the average of all features over a list of objects
     *  @param objects A list of objects (can be singlets or multiplets).
     *  @return An object vector where the first element is an object with the average features.
     */
    template<class T >
    static std::vector<T > average(const std::vector<T >& objects)
    {
        T object;
        std::vector<Matrix2D > featureMatrices = objects[0].features();
        for (int32 ind = 1; ind < objects.size(); ind ++) {
            for (int32 indF = 0; indF < objects[ind].features().size(); indF ++) {
                featureMatrices[indF] += objects[ind].features()[indF];

                if (ind == objects.size() - 1) {
                    object.add_feature(
                        objects[ind].names()[indF], 
                        featureMatrices[indF] / static_cast<MatrixElem >(objects.size()));
                }
            }
        }

        std::vector<T > vec;
        vec.push_back(object);

        return vec;
    };
};

}

#endif /* __AVERAGE_OBJECT_HXX__ */
