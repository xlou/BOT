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

#ifndef __TRACKING_PREDICTOR_HXX__
#define __TRACKING_PREDICTOR_HXX__

#include "FramePair.hxx"
#include "HammingLossFunction.hxx"
#include "CPLEXSolverSystem.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that predicts the tracking between a pair of frames using ILP
 *
 */
class TrackingPredictor
{
public:
	/*! Default constructor 
	 *
	 */
    TrackingPredictor()
    {
    };

	/*! Predict tracking
     *  @param framepair A FramePair object that contains the needed information
     *  @param weights The weights of features
     *  @param solution The return solution
     *  @param truth The ground truth, if not empty, this operation transforms to 
     *  solving the augmented problem for finding the most violated constraints.
	 */
    std::string operator()(
        const FramePair& framepair, 
        const std::vector<Matrix2D >& weights,
        std::vector<Matrix2D >& solution, 
        const std::vector<Matrix2D >& truth)
    {
        CPLEXSolverSystem solver;

        // is it the augmented optimization problem
        bool is_augmented = truth.size() != 0;

        // compute the losses, if necessary
        std::vector<Matrix2D > losses; 
        if (is_augmented) {
            HammingLossFunction lossFunc;
            lossFunc(truth, losses);
        }

        // formulate the coefficients in the objective function
        const std::vector<Event >& events = framepair.events();
        Matrix2D f;
        for (int32 ind = 0; ind < events.size(); ind ++) {
            Matrix2D tmp = events[ind].features() * weights[ind];
            if (is_augmented)
                tmp += losses[ind];
            f = (f.size() == 0) ? tmp : joinVertically(f, tmp);
        }

        // solve the ilp problem
        //std::cout << f << std::endl;
        const Matrix2D& Aeq = framepair.constraint();
        Matrix2D beq = Matrix2D(Shape2D(Aeq.shape(0), 1), static_cast<MatrixElem >(1));
        Matrix2D MatrixNull;
        Matrix2D x;
        std::string msg = solver.solve_bilp(
            -f,                         /* the objective */
            MatrixNull, MatrixNull,     /* the inequality constraints */
            Aeq, beq,                   /* the equality constraints */
            MatrixNull,                 /* the initial solution */
            x);                         /* the solution */

        // fold ilp solution
        solution = framepair.fold_solution(x);

        return msg;
    };
};

}

#endif /* __TRACKING_PREDICTOR_HXX__ */
