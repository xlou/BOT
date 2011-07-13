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

#ifndef __HAMMING_LOSS_FUNCTION_HXX__
#define __HAMMING_LOSS_FUNCTION_HXX__

#include <iostream>
#include <functional>
#include "TypeDefinition.hxx"
#include "vigra/matrix.hxx"
#include "VigraSTLInterface.hxx"

using namespace vigra::linalg;

namespace bot {

/*! A class which computes the loss or the coefficient of the loss function 
 *  as a penalty function of a prediction.
 *
 */
class HammingLossFunction
{
public:
	/*! Default constructor 
	 *
	 */
    HammingLossFunction() {};

	/*! Computes the normalized Hamming distance between two binary vectors.
     *  @param z The ground truth as a vector of Matrix2D
     *  @param z_hat The prediction as a vector of Matrix2D
     *  @return MatrixElem The loss as a MatrixElem value
	 */
    MatrixElem loss(
        const Solution& z, 
        const Solution& z_hat) const
    {

		MatrixElem l;
		operator()(z, z_hat, l);

        return l;
    };

	/*! Computes the normalized Hamming distance between two binary vectors.
     *  @param z The ground truth as a vector of Matrix2D
     *  @param z_hat The prediction as a vector of Matrix2D
     *  @param loss The return loss as a MatrixElem value
	 */
    void operator()(
        const Solution& z, 
        const Solution& z_hat,
        MatrixElem& l) const
    {
        l = 0;
        MatrixElem sum = 0;
        for (int32 ind = 0; ind < z.size(); ind ++) {
            Matrix2D tmp = abs(z_hat[ind] - z[ind]);
            l += std::accumulate(
                tmp.begin(), 
                tmp.end(), 
                static_cast<MatrixElem >(0));
            sum += z[ind].size(0);
        }

        l /= sum;
    };

	/*! Computes the coefficients of the loss as a penalty function of a 
     *  prediction. The loss, given a prediction z_hat, is [numOnes/sum + <b, z_hat>].
     *  @param z The ground truth as a vector of Matrix2D
     *  @pram b The coefficients of the loss as a matrix
	 */
    void operator()(
        const std::vector<Matrix2D >& z, 
        std::vector<Matrix2D >& b) const
    {
        MatrixElem numOnes = VigraSTLInterface::accumulate(z);

        MatrixElem sum = 0;
        for (int32 ind = 0; ind < z.size(); ind ++)
            sum += z[ind].size(0);

        b.clear();

        for (int32 ind = 0; ind < z.size(); ind ++) {
			Matrix2D one(Shape2D(z[ind].size(0), 1), static_cast<MatrixElem>(1));
            b.push_back((one - z[ind]*static_cast<MatrixElem>(2))/sum);
		}
    };
};

}

#endif /* __HAMMING_LOSS_FUNCTION_HXX__ */
