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

#ifndef __TYPE_DEFINITION_HXX__
#define __TYPE_DEFINITION_HXX__

#include <vector>
#include <map>
#include "vigra/multi_array.hxx"

namespace bot {

typedef int int32;
typedef unsigned short uint16;

typedef double MatrixElem;
typedef vigra::MultiArray<2, MatrixElem > Matrix2D;
typedef Matrix2D::difference_type Shape2D;
typedef Matrix2D::view_type View2D;

#define BASE_ID             1

#define LARGEST_DISTANCE    100000
#define DEFAULT_NEIGHBORS   3

typedef std::pair<int32, int32 > IDPair;
typedef std::pair<int32, int32 > LabelPair;

typedef std::pair<std::string, std::string > FeatureMeasure;
typedef std::vector<FeatureMeasure > MeasureSetting;

typedef std::pair<int32, int32 > EventPairing;

typedef std::vector<std::pair<int32, int32 > > Hypotheses;

typedef std::pair<Matrix2D, Matrix2D > Mapping;

typedef std::vector<Matrix2D > Solution;

struct LabelAssociation {
    std::string name;
    Matrix2D source;
    Matrix2D target;
};

typedef std::vector<LabelAssociation > LabelAssociations;

}
#endif /* __TYPE_DEFINITION_HXX__ */
