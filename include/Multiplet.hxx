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

#ifndef __MULTIPLET_HXX__
#define __MULTIPLET_HXX__

#include "Object.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that represents a multiplet, i.e. a combination of single 
 *  segments (aka singlets)
 *
 */
class Multiplet : public Object
{
public:
	/*! Default constructor 
	 *
	 */
    Multiplet()
    {
    };

	/*! Constructor with IDs, pixels, labels and values specified
     *  @param id The id of this multiplet
     *  @param comp_id_1 The id of the first component
     *  @param comp_id_2 The id of the second component
     *  @param pixels The pixels (coordinates) that this multiplet contains
     *  @param values The values (intensities) that this multiplet contains
     *  @param labels The segmentation labels that this multiplet contains
	 */
    Multiplet(
        const int32 id, const int32 comp_id_1, const int32 comp_id_2, 
        const Matrix2D& pixels, const Matrix2D& values, const Matrix2D& labels)
        : Object(id, pixels, values, labels)
    {
        components_.push_back(std::min(comp_id_1, comp_id_2));
        components_.push_back(std::max(comp_id_1, comp_id_2));
    };
};

typedef std::vector<Multiplet > Multiplets;
typedef std::vector<Multiplets > MultipletsSequence;

}

#endif /* __MULTIPLET_HXX__ */
