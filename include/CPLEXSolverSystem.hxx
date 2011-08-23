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

#ifndef __CPLEX_SOLVER_SYSTEM_HXX__
#define __CPLEX_SOLVER_SYSTEM_HXX__

#include "TypeDefinition.hxx"
#include "SolverSystem.hxx"
#include "vigra/matrix.hxx"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using namespace vigra::linalg;

namespace bot 
{

/*! A class that solves ILP/QP problems using CPLEX
 *
 */
class CPLEXSolverSystem : public SolverSystem 
{
public:
    /*! Default constructor, initialize the concert environment
     *  
     */
    CPLEXSolverSystem()
    {
        env.setOut(env.getNullStream());
    };

    /*! Destructor, end the concert environment
     *  
     */
    ~CPLEXSolverSystem()
    {
        env.end();
    };

    /*! Setup the binary ilp objective function
     *  @param f The coefficients of the objective function
     *  @param vars The random variables
     *  @param model The model where the optimization problem will be added
     */
    void setup_lp_objective(
        const Matrix2D& f, 
        const IloNumVarArray& vars,
        IloModel& model) const
    {
        IloObjective objective(env);
        for (int32 elem = 0; elem < f.size(); elem ++) 
            objective.setLinearCoef(vars[elem], static_cast<IloNum >(f[elem]));
        model.add(objective);
    };

    /*! Extract the solution
     *  @param cplex The CPLEX solver object
     *  @param vars The random variables
     *  @return x The solution as a matrix
     */
    void extract_solution(
        const IloCplex& cplex, 
        const IloNumVarArray& vars,
        Matrix2D& x) const
    {
        IloNumArray vals(env);
		std::cout << "getting values..." << std::endl;
        cplex.getValues(vals, vars);
        x.reshape(Shape2D(vals.getSize(), 1), static_cast<MatrixElem >(0));
		std::cout << "extracting " << vals.getSize() << " variables" << std::endl;
        for (int32 ind = 0; ind < vals.getSize(); ind ++)
            x[ind] = static_cast<MatrixElem >(vals[ind]);
    };

    /*! Setup linear constraints, i.e.
     *      lb <= Ax <= ub
     *  @param A The constraint matrix
     *  @param lb The lower bound
     *  @param up The upper bound
     *  @param vars The random variables
     *  @param model The cplex optimziation model
     */
    void setup_linear_constraints(
        const Matrix2D& A, const Matrix2D& lb, const Matrix2D& ub, 
        const IloNumVarArray& vars, IloModel& model) const
    {
        if (A.size() == 0)
            return ;

        IloRangeArray constraints(env);

        // add constraints one by one
        for (int32 row = 0; row < A.shape(0); row ++) {
            // set the coefficients
            IloRange constraint(env, -IloInfinity, IloInfinity);

            for (int32 col = 0; col < A.shape(1); col ++) {
                if (A(row, col) == 0)
                    continue;

                constraint.setLinearCoef(vars[col], static_cast<IloNum >(A(row, col)));
            }

            // add bounds
            if (lb.shape(0) == A.shape(0))
                constraint.setLB(lb(row, 0));
            if (ub.shape(0) == A.shape(0))
                constraint.setUB(ub(row, 0));

            // add to the constraint array
            constraints.add(constraint);
        }

        model.add(constraints);
    };

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
    std::string solve_bilp(
        const Matrix2D& f, 
        const Matrix2D& Aineq, const Matrix2D& bineq, 
        const Matrix2D& Aeq, const Matrix2D& beq, 
        const Matrix2D& x0, 
        Matrix2D& x) const
    {
        std::string msg("Error occured during cplex problem formulation");
        try {
            // create the model and the random variables
            IloModel model(env);
            IloNumVarArray vars(env, f.size(), 0, 1, ILOINT);

			std::cout << "setting up lp: " << f.size(0) << " variables, "
			          << bineq.size(0) << " inequalities, "
			          << beq.size(0) << " equalities" << std::endl;

            // add the constraints
            Matrix2D lb;
            setup_linear_constraints(Aineq, lb, bineq, vars, model);   // inequality
            setup_linear_constraints(Aeq, beq, beq, vars, model);           // equality

            // add the objective functions
            setup_lp_objective(f, vars, model);

			// make sure all variables have been added to the model
			for (int i = 0; i < f.size(0); i++)
				model.add(IloRange(env, vars[i].getLB(), vars[i], vars[i].getUB()));

			std::cout << "solving integer linera program" << std::endl;

            // call cplex to solve it
            IloCplex cplex(model);
            cplex.solve();

			std::cout << "done: " << cplex.getStatus() << std::endl;
			std::cout << "querying result status" << std::endl;

            // get cplex status
            if (cplex.getStatus() == IloAlgorithm::Optimal)
                msg = "Optimal solution obtained";
            else {
                msg = "*Warning* Fail to obtain the optimal solution";
                std::cout << "cplex.getStatus() gives: " << cplex.getStatus() << std::endl;
            }

			std::cout << "extracting solution" << std::endl;

            // extract the solution
            extract_solution(cplex, vars, x);
        }
        catch (IloException& e) {
            msg = std::string("*Warning* Concert exception caught: ") + e.getMessage();
        }
        catch (...) {
            msg = "*Warning* Unknown exception caught";
        }

        return msg;
    };

    /*! Setup the quadratic programming objective function
     *  @param H The coefficients of the quadratic objective
     *  @param f The coefficients of the linear objective
     *  @param vars The random variables
     *  @param model The model where the optimization problem will be added
     */
    void setup_qp_objective(
        const Matrix2D& H, const Matrix2D& f, 
        const IloNumVarArray& vars,
        IloModel& model) const
    {
        IloNumExpr expr(env);

        // add the quadratic terms
        for (int32 row = 0; row < H.shape(0); row ++) {
            for (int32 col = row; col < H.shape(1); col ++) {
                MatrixElem h = H(row, col);
                if (h == 0)
                    continue ;

                if (col == row) 
                    expr = expr + static_cast<IloNum >(0.5*h)*vars[row]*vars[row];
                else
                    expr = expr + static_cast<IloNum >(h)*vars[row]*vars[col];
            }
        }

        // add the linear terms
        for (int32 elem = 0; elem < f.size(); elem ++) {
            if (f[elem] == 0)
                continue;

            expr = expr + static_cast<IloNum >(f[elem])*vars[elem];
        }

        // add to the model
        model.add(IloObjective(env, expr));
    };

	/*! Setup the bounds
     *  @param vars The array of random variables
     *  @param lb A Matrix2D object that contains the lower bounds
     *  @param ub A Matrix2D object that contains the upper bounds
     */
    void setup_bounds(
        IloNumVarArray& vars,
        const Matrix2D& lb, const Matrix2D& ub) const
    {
        // set the lower bounds
        IloNumArray lbAry(env, vars.getSize());
        for (int32 ind = 0; ind < lbAry.getSize(); ind ++) 
            lbAry[ind] = lb.size() == 0 ? -IloInfinity : lb[ind];

        // set the upper bounds
        IloNumArray ubAry(env, vars.getSize());
        for (int32 ind = 0; ind < ubAry.getSize(); ind ++) 
            ubAry[ind] = ub.size() == 0 ? IloInfinity : ub[ind];

        vars.setBounds(lbAry, ubAry);
    };

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
    std::string solve_qp(
        const Matrix2D& H, const Matrix2D& f, 
        const Matrix2D& Aineq, const Matrix2D& bineq, 
        const Matrix2D& Aeq, const Matrix2D& beq, 
        const Matrix2D& lb, const Matrix2D& ub, 
        const Matrix2D& x0, Matrix2D& x) const
    {
        std::string msg("Error occured during cplex problem formulation");
        try {
            // create the model and the random variables
            IloModel model(env);
            IloNumVarArray vars(env, f.size(), -IloInfinity, IloInfinity, ILOFLOAT);
            setup_bounds(vars, lb, ub);

            // add the constraints
            Matrix2D empty_matrix;
            setup_linear_constraints(Aineq, empty_matrix, bineq, vars, model);   // inequality
            setup_linear_constraints(Aeq, beq, beq, vars, model);           // equality

            // add the objective functions
            setup_qp_objective(H, f, vars, model);

			// make sure all variables have been added to the model
			for (int i = 0; i < f.size(0); i++)
				model.add(IloRange(env, vars[i].getLB(), vars[i], vars[i].getUB()));

            // call cplex to solve it
            IloCplex cplex(model);
            cplex.solve();

            if (cplex.getStatus() == IloAlgorithm::Optimal)
                msg = "Optimal solution obtained";
            else {
                msg = "*Warning* Fail to obtain the optimal solution";
                std::cout << "cplex.getStatus() gives: " << cplex.getStatus() << std::endl;
            }
            extract_solution(cplex, vars, x);
        }
        catch (IloException& e) {
            msg = std::string("*Warning* Concert exception caught: ") + e.getMessage();
        }
        catch (...) {
            msg = "*Warning* Unknown exception caught";
        }

        return msg;
    };

private:
    IloEnv env;
};

}

#endif /* __CPLEX_SOLVER_SYSTEM_HXX__ */
