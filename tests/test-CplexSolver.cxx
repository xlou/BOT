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

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
int
main (int argc, char **argv)
{
    IloEnv env;
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        vars.add(IloNumVar(env, 0.0, 40.0));
        vars.add(IloNumVar(env));
        vars.add(IloNumVar(env));
        model.add(IloMaximize(env, vars[0] + 2 * vars[1] + 3 * vars[2]));
        model.add( - vars[0] + vars[1] + vars[2] <= 20);
        model.add( vars[0] - 3 * vars[1] + vars[2] <= 30);
        IloCplex cplex(model);
        if ( !cplex.solve() ) {
            env.error() << "Failed to optimize LP." << endl;
            throw(-1);
        }
        IloNumArray vals(env);
        env.out() << "Solution status = " << cplex.getStatus() << endl;
        env.out() << "Solution value = " << cplex.getObjValue() << endl;
        cplex.getValues(vals, vars);
        env.out() << "Values = " << vals << endl;
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();

    return 0;
}