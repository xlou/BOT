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

#include "../include/NearestNeighborGenerator.hxx"
#include <vector>

using namespace bot;

int main()
{
    MatrixElem array1[] = {81,91,13,91,63,10,28,55,96,96,16,97,96,49,80,14,42,92,79,96};
    MatrixElem array2[] = {44,38,77,80,19,49,45,65,71,75,28,68,66,16,12,50};

    Matrix2D pts1(Matrix2D::difference_type(10, 2), 0.0);
    Matrix2D pts2(Matrix2D::difference_type(8, 2), 0.0);
    for (int32 ind = 0; ind < pts1.size(); ind ++) 
        pts1[ind] = array1[ind];
    for (int32 ind = 0; ind < pts2.size(); ind ++) 
        pts2[ind] = array2[ind];

    // test symmetrical neighborhood search
    NearestNeighborGenerator nnGenSym(3, 40, true);
    std::cout << "**** test symmetrical neighborhood search ****" << std::endl;
    std::cout << nnGenSym(pts1, pts1) << std::endl;

    // test asymmetrical neighborhood search
    NearestNeighborGenerator nnGenAsym(3, 40, false);
    std::cout << "**** test asymmetrical neighborhood search ****" << std::endl;
    std::cout << nnGenAsym(pts1, pts2) << std::endl;

    return 0;
}
