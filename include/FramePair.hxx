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

#ifndef __FRAME_PAIR_HXX__
#define __FRAME_PAIR_HXX__

#include <vector>
#include <vigra/matrix.hxx>
#include "TypeDefinition.hxx"
#include "Event.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! This class contains all information for a pair of neighboring frames, 
 *  such as the hypotheses, the feature matrix and the coupling constraint.
 *
 */
class FramePair {
public:
    /*! Default constructor
     *
     */
    FramePair() {};

    /*! Add and event, including hypotheses and feature matrix
     *  @param name The name of the event
     *  @param pairing The pairing of the event
     *  @param hypotheses Pairwise hypotheses
     *  @param features The feature matrix for these hypotheses
     *  @param constraint The coupling constraints
     */
    void add_event(
        const std::string& name, 
        const EventPairing& pairing, 
        const Hypotheses& hypotheses,
        const Matrix2D& features, 
        const Matrix2D& constraint)
    {
        events_.push_back(Event(name, pairing, hypotheses, features));
        constraint_ = (constraint_.size() == 0) ? 
            constraint : joinHorizontally(constraint_, constraint);
    };

    /*! Return a constant reference to the events
     *  @return A std::vector<Event > object as the reference
     */
    const std::vector<Event >& events() const
    {
        return events_;
    };

    /*! Return a reference to the events
     *  @return A std::vector<Event > object as the reference
     */
    std::vector<Event >& events()
    {
        return events_;
    };
    
    /*! Return a constant reference to the constraint matrix
     *  @return A Matrix2D object as the constraint matrix
     */
    const Matrix2D& constraint() const
    {
        return constraint_;
    };

    /*! Return a constant reference to the solution
     *  @return A constant reference to the solution
     */
    const Solution& solution() const
    {
        return solution_;
    };

    /*! Return a reference to the solution
     *  @return A reference to the solution
     */
    Solution& solution() 
    {
        return solution_;
    };

    /*! Fold the ilp solution (a single matrix) into a vector of matrices
     *  @param x The ilp solution
     *  @return A std::vector<Matrix2D > object as the folded solution
     */
    const std::vector<Matrix2D > fold_solution(const Matrix2D& x) const
    {
        std::vector<Matrix2D > solution;
        const std::vector<Event >& vec = events();
        int32 sum = 0;
        for (int32 ind = 0; ind < vec.size(); ind ++) {
            Shape2D p(sum, 0), q(sum+vec[ind].size(), 1);
            solution.push_back(Matrix2D(x.subarray(p, q)));
            sum += vec[ind].size();
        }

        return solution;
    };

    /*! Unfold the solution (a vector of matrices) into a single matrix
     *  @param solution The folded solution
     *  @return A Matrix2D object as the unfolded solution
     */
    Matrix2D unfold_solution(const std::vector<Matrix2D >& solution) 
    {
        Matrix2D x;
        for (int32 ind = 0; ind < solution.size(); ind ++) 
            x = joinVertically(x, solution[ind]);

        return x;
    };

    /*! Return an empty solution
     *  @return A std::vector<Matrix2D > object as the empty solution
     */
    std::vector<Matrix2D > empty_solution() 
    {
        std::vector<Matrix2D > solution;
        for (int32 ind = 0; ind < events_.size(); ind ++) {
            Matrix2D mat(
                Matrix2D::difference_type(events_[ind].hypotheses().size(), 1), 
                static_cast<MatrixElem >(0));
            solution.push_back(mat);
        }

        return solution;
    };

    /*! Compute the compatibility vector from a solution
     *  @param x A Solution object as the solution
     *  @return A Matrix2D object as the compatibility vector
     */
    Matrix2D compatibility_vector(const Solution& z) const
    {
        Matrix2D mat;
        for (int32 ind = 0; ind < events().size(); ind ++) {
            mat = (mat.size() == 0) ? 
                transpose(events()[ind].features()) * z[ind] : 
                joinVertically(mat, transpose(events()[ind].features()) * z[ind]);
        }

        return mat;
    };
protected:
    std::vector<Event > events_;
    Matrix2D constraint_;
    Solution solution_;
};

}

#endif /* __FRAME_PAIR_HXX__ */
