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

#ifndef __HYPOTHESIS_GENERATOR_HXX__
#define __HYPOTHESIS_GENERATOR_HXX__

#include "Object.hxx"
#include "VigraSTLInterface.hxx"
#include "NearestNeighborGenerator.hxx"

namespace bot 
{

/*! This struct compares two std::pair type objects first by the
 *  second element then by the first one.
 */
template<class T >
struct PairLessThan 
{
    bool operator()(const std::pair<T, T > &left, const std::pair<T, T > &right) 
    {
        if (left.second == right.second)
            return left.first < right.first;

        return left.second < right.second;
    }
};

/*! A class that generates hypotheses of pairwise assignments.
 *
 */
class HypothesisGenerator : public NearestNeighborGenerator
{
public:
	/*! Constructor with neighborhood number and maximal distance specified
	 *
	 */
    HypothesisGenerator(const Index k, const CoordElem d_max)
    {
        k_ = k;
        d_max_ = d_max;

    };

	/*! Overridden operator(), generate hypytheses between two lists of objects. 
     *  Each list can be a list of singlets or multiplets.
	 *  @param objects1 A vector of objects as the source
     *  @param objects2 A vector of objects as the target
     *  @param full If true, create a fully connected bipartite graph
     *  @return A vector of hypothesis (pairs of indices of objects)
	 */
    template<class T1, class T2 >
    Hypotheses operator()(
        const std::vector<T1 >& objects1, 
        const std::vector<T2 >& objects2,
        const bool full = false) 
    {
        Hypotheses hypotheses;
//	if (objects1.size() == 0 || objects2.size() == 0)
//		return hypotheses;

        if (full) {
            for (Index ind1 = 0; ind1 < objects1.size(); ind1 ++)
                for (Index ind2 = 0; ind2 < objects2.size(); ind2 ++)
                    hypotheses.push_back(std::make_pair(ind1, ind2));
        }
        else {
            // search for id pairs in a neighborhood (kNN + distance thresholding)
            CoordMatrix centers1 = VigraSTLInterface::get_centers<T1 >(objects1);
            CoordMatrix centers2 = VigraSTLInterface::get_centers<T2 >(objects2);
            hypotheses = NearestNeighborGenerator::operator()(centers1, centers2);
        }

        // sort them first by the right element and then by the left element
        std::sort(hypotheses.begin(), hypotheses.end(), PairLessThan<int32 >());

        return hypotheses;
    };
};

}

#endif /* __HYPOTHESIS_GENERATOR_HXX__ */
