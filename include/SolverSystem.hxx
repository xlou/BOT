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

#ifndef __SOLVER_SYSTEM_HXX__
#define __SOLVER_SYSTEM_HXX__

#include <iostream>
#include "TypeDefinition.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! The base class for ILP/QP solver
 *
 */
class SolverSystem
{
public:
    /*! Solve a binary integer linear programminig (binary ilp) problem with 
     *  given equality and inequality constraints. In particular, it solves 
     *  a problem as follows:
     *          max     f'*x
     *          s.t.    Aineq*x <=  bineq
     *                  Aeq*x   =   beq
     *                  each variable in x is binary
     *  @param f The coefficient of the objective function
     *  @param Aineq The inequality constraints: Aineq * x <= bineq
     *  @param bineq The inequality constraints: Aineq * x <= bineq
     *  @param Aeq The equality constraints: Aineq * x = bineq
     *  @param beq The equality constraints: Aineq * x = bineq
     *  @param x0 The initial solution
     *  @param x The solution
     *  @param msg The return message
     */
    virtual std::string solve_bilp(
        const Matrix2D& f, 
        const Matrix2D& Aineq, const Matrix2D& bineq, 
        const Matrix2D& Aeq, const Matrix2D& beq, 
        const Matrix2D& x0, Matrix2D& x) const = 0;

	/*! Solve a quadaratic programminig (qp) problem with given equality and 
     *  inequality constraints and bounds. In particular, it solve:
     *      min     0.5*x'*H*x+f*x
     *      st.     Aineq*x <= bineq
     *              Aeq*x    = beq
     *              lb <= x <= ub
     *  @param H Double matrix for objective function (quadratic term)
     *  @param f Double matrix (vector) for objective function (linear term)
     *  @param Aineq The inequality constraints: Aineq * x <= bineq
     *  @param bineq The inequality constraints: Aineq * x <= bineq
     *  @param Aeq The equality constraints: Aineq * x = bineq
     *  @param beq The equality constraints: Aineq * x = bineq
     *  @param lb The lower bound
     *  @param ub The upper bound
     *  @param x0 The initial solution
     *  @param x The solution
     *  @return A std::string as the message of the solution status
	 */
    virtual std::string solve_qp(
        const Matrix2D& H, const Matrix2D& f, 
        const Matrix2D& Aineq, const Matrix2D& bineq, 
        const Matrix2D& Aeq, const Matrix2D& beq, 
        const Matrix2D& lb, const Matrix2D& ub, 
        const Matrix2D& x0, Matrix2D& x) const = 0;
};

}

#endif /* __SOLVER_SYSTEM_HXX__ */
