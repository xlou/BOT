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

#ifndef __TRACKING_TRAINER_HXX__
#define __TRACKING_TRAINER_HXX__

#include "VigraSTLInterface.hxx"
#include "FramePair.hxx"
#include "HammingLossFunction.hxx"
#include "CPLEXSolverSystem.hxx"
#include "TrackingPredictor.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that learns the optimal weights from true tracking results
 *
 */
class TrackingTrainer
{
public:
	/*! Default constructor with default initiliazation of the parameters
	 *
	 */
    TrackingTrainer()
    {
        epsilon_ = 2.0e-4;
        eta_ = 1e-10;
        lambda_ = 1e-8;
        max_iter_ = 500; 
    };

	/*! Constructor with specified parameters
	 *
	 */
    TrackingTrainer(
        const double epsilon, const double eta,
        const double lambda, const double max_iter) :
    epsilon_(epsilon), eta_(eta), lambda_(lambda), max_iter_(max_iter)
    {
    };

    /*! Fold the weights into a vector of matrices
     *  @param W The weight vector
     *  @return A std::vector<Matrix2D > object as the folded weights
     */
    void fold_weights(const Matrix2D& W, std::vector<Matrix2D >& weights) 
    {
        int32 sum = 0;
        for (int32 ind = 0; ind < weights.size(); ind ++) {
            Shape2D p(sum, 0), q(sum+weights[ind].size(), 1);
            sum += weights[ind].size();
            weights[ind] = Matrix2D(W.subarray(p, q));
        }
    };

    /*! Update the weights by solving the QP problem
     *  @param W The weight vector
     *  @return A std::vector<Matrix2D > object as the updated weights
     */
    void update_weights(
        const SolverSystem& solver, 
        const Matrix2D& matA, const Matrix2D& matB, 
        Matrix2D& W)
    {
        Matrix2D Alpha, empty_matrix;
        Matrix2D Aeq(Shape2D(1, matB.size()), static_cast<MatrixElem >(1));
        Matrix2D beq(Shape2D(1, 1), static_cast<MatrixElem >(1));
        Matrix2D lb(Shape2D(matB.size(), 1), static_cast<MatrixElem >(0));

        std::string message = solver.solve_qp(
            transpose(matA) * matA / lambda_, -matB, 
            empty_matrix, empty_matrix, 
            Aeq, beq, 
            lb, empty_matrix, 
            empty_matrix, Alpha);

		if (verbose)
			std::cout << "Solver returned: " << message << std::endl;

        W = -matA * Alpha / lambda_;
    };

	/*! Learn the weights
     *  @param training A TrainingData object that contains the training data
     *  @param framepairs All the frame pairs
     *  @param weights The initial weights and the return one
     *  @param verbose If true, output information of the intermediate steps 
     *  @return A std::string object as the return message
	 */
    std::string operator()(
        const TrainingData& training, 
        const std::vector<FramePair >& framepairs, 
        std::vector<Matrix2D >& weights, 
        const bool verbose = false)
    {
        std::string msg;
        CPLEXSolverSystem solver;

        // the loss and the inference method
        HammingLossFunction lossFunc;
        TrackingPredictor predictor;

        // some intermediate results
        Matrix2D matA, matB, matW;      //       
        Matrix2D matLosses;             // intermediate losses 
        Matrix2D matEpsilons;             // intermediate approximate gaps 
        Matrix2D J;                // original object function
        double J_previous = 1e10;            // last approximation
        
        // enter the iterations
        const std::vector<int32 >& times = training.times();
        int32 nTr = times.size();
        const std::vector<Solution >& solutionsTr = training.solutions();
        int32 iter = 0;
        const std::vector<Matrix2D > null_vector;
        while (iter < max_iter_) {
            if (verbose)
                std::cout << "****Training iteration #" << iter << "****" << std::endl;

            // determine most violated constraints
            std::vector<Solution > solutionsPred;
            for (int32 indTr = 0; indTr < nTr; indTr ++) {
                const FramePair& framepair = framepairs[times[indTr]];
                std::vector<Matrix2D > solution;
                predictor(framepair, weights, solution, solutionsTr[indTr]);
                solutionsPred.push_back(solution);
            }
            if (verbose)
                std::cout << "\t\tMost violated constraints determined." << std::endl;

            // compute the gradients and offsets
            Matrix2D A, B;
            Matrix2D W = VigraSTLInterface::unfold(weights);
            for (int32 indTr = 0; indTr < nTr; indTr ++) {
                const FramePair& framepair = framepairs[times[indTr]];
                Matrix2D PhiTr = framepair.compatibility_vector(solutionsTr[indTr]);
                Matrix2D Phi = framepair.compatibility_vector(solutionsPred[indTr]);
                Matrix2D Psi = PhiTr - Phi;
                if (A.size() == 0) {
                    A = -Psi;
                    B.reshape(Shape2D(1, 1), static_cast<MatrixElem >(0));
                    B[0] = lossFunc.loss(solutionsTr[indTr], solutionsPred[indTr]) - dot(W, Psi);
                }
                else {
                    A -= Psi;
                    B[0] += lossFunc.loss(solutionsTr[indTr], solutionsPred[indTr]) - dot(W, Psi);
                }
            }
            A /= static_cast<MatrixElem >(nTr);
            B /= static_cast<MatrixElem >(nTr);
			B -= dot(W, A);
            if (verbose)
                std::cout << "\t\tGradients and offsets computed." << std::endl;

            // solve the dual problem with quadratic programming
            matA = (matA.size() == 0) ? A : joinHorizontally(matA, A);
            matB = (matB.size() == 0) ? B : joinHorizontally(matB, B);
            update_weights(solver, matA, matB, W);
            matW = (matW.size() == 0) ? W : joinHorizontally(matW, W);
            if (verbose)
                std::cout << "\t\tWeights updated (L2 regularizer, QP problem solved)." << std::endl;

            // perform prediction using the new weights
            fold_weights(W, weights);
            Matrix2D Losses(Shape2D(nTr, 1), static_cast<MatrixElem >(0));
			// do this for every training sample to get the current value of
			// l(x,y,w_t) -- this will be used later to get the current value of
			// J(w_t)
            for (int32 indTr = 0; indTr < nTr; indTr ++) {

                const FramePair& framepair = framepairs[times[indTr]];

				// get y* = argmax_y' l(x_i,y_i,w_t)
                std::vector<Matrix2D > solution;
                predictor(framepair, weights, solution, solutionsTr[indTr]);

				// compute l_i = l(x_i,y*,w_t)
                Losses[indTr] = lossFunc.loss(solutionsTr[indTr], solution);

				Matrix2D PhiTr = framepair.compatibility_vector(solutionsTr[indTr]);
				Matrix2D Phi   = framepair.compatibility_vector(solution);
				Matrix2D Psi   = Phi - PhiTr;
				MatrixElem diff = dot(W,Psi);

				Losses[indTr] += diff;

				std::cout << "*****************************************" << std::endl;
				//std::cout << "W" << std::endl << W << std::endl;
				//std::cout << "PhiTr" << std::endl << PhiTr << std::endl;
				//std::cout << "Phi" << std::endl << Phi << std::endl;
				//std::cout << "Psi" << std::endl << Psi << std::endl;
				//std::cout << "diff" << std::endl << (transpose(W)*Psi) << std::endl;
				std::cout << "loss" << std::endl
				          << diff << " + "
				          << lossFunc.loss(solutionsTr[indTr], solution) << " = "
						  << Losses[indTr] << std::endl;
				std::cout << "*****************************************" << std::endl;
            }
            matLosses = (matLosses.size() == 0) ? Losses : joinHorizontally(matLosses, Losses);
            if (verbose)
                std::cout << "\t\tEmpirical loss updated." << std::endl;

            // compute the approximation gap
            double reg = squaredNorm(W)/2.0;
            Matrix2D P = transpose(matA) * W + transpose(matB);
            double approxJ = lambda_ * reg + *std::max_element(P.begin(), P.end());
            double J_ = lambda_ * reg + std::accumulate(Losses.begin(), Losses.end(), static_cast<MatrixElem >(0)) / static_cast<MatrixElem >(Losses.size());
            J = (J.size() == 0) ? 
                Matrix2D(Shape2D(1, 1), static_cast<MatrixElem >(J_)) : 
                joinHorizontally(J, Matrix2D(Shape2D(1, 1), static_cast<MatrixElem >(J_)));
            double gap = *std::min_element(J.begin(), J.end()) - approxJ;
            matEpsilons = (matEpsilons.size() == 0) ? 
                Matrix2D(Shape2D(1, 1), static_cast<MatrixElem >(gap)) : 
                joinHorizontally(matEpsilons, Matrix2D(Shape2D(1, 1), static_cast<MatrixElem >(gap)));
            if (verbose)
                std::cout << "\t\tApproximate gap computed." << std::endl;

            // determine whether to stop
            if (gap < 0) 
                std::cerr << "*Warning* Possible lower bound error: approximation gap = " << gap << std::endl;

            if (gap < epsilon_) {
                msg = "Precision reaches the given epsilon";
                if (verbose)
                    std::cout << "\t\tPrecision reaches the given epsilon." << std::endl;
                break ;
            }

            if (std::abs(J_ - J_previous) < eta_) {
                msg = "Approximation improvement reaches the given eta";
                if (verbose)
                    std::cout << "\t\tApproximation improvement reaches the given eta." << std::endl;
                break ;
            }
            
            // next iteration
            J_previous = J_;
            iter ++;
        }

        // store the intermediate results
        matW_ = matW;
        matEpsilons_ = matEpsilons;
        matLosses_ = matLosses;
        iter_ = iter;
        if (iter == max_iter_)
            msg = "Maximum iteration reached";

        return msg;
    };

    /*! Get a constant reference to the weights vector
     *  @return A constant reference to the weights vector
     */
    const Matrix2D& weights() const
    {
        return matW_;
    };
    
    /*! Get a constant reference to the epsilon vector
     *  @return A constant reference to the epsilon vector
     */
    const Matrix2D& epsilons() const
    {
        return matEpsilons_;
    };
    
    /*! Get a constant reference to the losses vector
     *  @return A constant reference to the losses vector
     */
    const Matrix2D& losses() const
    {
        return matLosses_;
    };

    /*! Return the number of iterations before convergence
     *  @return A int32 value as the number of iterations
     */
    const int32 iter() const
    {
        return iter_;
    };
protected:
    double epsilon_;
    double eta_;
    double lambda_;
    int32 max_iter_;
    int32 iter_;
    Matrix2D matW_;
    Matrix2D matEpsilons_;
    Matrix2D matLosses_;
};

}

#endif /* __TRACKING_TRAINER_HXX__ */
