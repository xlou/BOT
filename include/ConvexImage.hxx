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

#ifndef __CONVEX_IMAGE_HXX__
#define __CONVEX_IMAGE_HXX__

#include "vigra/matrix.hxx"
#include "vigra/polygon.hxx"
#include "TypeDefinition.hxx"
#include "VigraSTLInterface.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! The class computes the covnex hull a given list of points and outpus another list of points that locate inside the convex hull.
 *
 */
class ConvexImage 
{
public:
    /*! Return a binary mask image covering the convex hull of a given point cloud
     *  @return A binary mask image covering the convex hull
     */
    static Matrix2D convex_hull(const Matrix2D& pts)
    {
        //typedef vigra::ArrayVector<vigra::Point2D > t_pt_array;
        typedef std::vector<Matrix2D > t_pt_array;

        // create array
        t_pt_array pt_ary;
        for (int i = 0; i < pts.shape(0); i ++) {
            Matrix2D pt(Matrix2D::difference_type(1, 2), 0.0);
            pt[0] = pts(i, 0);
            pt[1] = pts(i, 1);
            pt_ary.push_back(pt);
        }

        // compute the convex hull
        t_pt_array hull;
        vigra::convexHull<t_pt_array, t_pt_array >(pt_ary, hull);

        // convert to Matrix2D
        Matrix2D convex_hull = 
            VigraSTLInterface::vector_to_matrix<Matrix2D, MatrixElem >(hull);

        return convex_hull;
    };

    /*! Return if a point is in the convex hull
     *  @return A bool value, true when the point is in the convex hull
     */
    static bool inside(const MatrixElem x, const MatrixElem y, const Matrix2D& convex_hull)
    {
        Matrix2D pt(Matrix2D::difference_type(1, 2), 0.0);
        pt[0] = x;
        pt[1] = y;
        return inside(pt, convex_hull);
    }

    /*! Return if a point is in the convex hull
     *  @return A bool value, true when the point is in the convex hull
     */
    static bool inside(const Matrix2D& pt, const Matrix2D& convex_hull)
    {
        Matrix2D v0(rowVector(convex_hull, convex_hull.shape(0)-1) - rowVector(convex_hull, convex_hull.shape(0)-2));
        v0 = v0 / norm(v0);
        
        for (int32 ind = 0; ind < convex_hull.shape(0)-1; ind ++) {
            Matrix2D v1(rowVector(convex_hull, ind+1) - rowVector(convex_hull, ind));
            v1 = v1 / norm(v1);

            Matrix2D v2(pt - rowVector(convex_hull, ind));
            if (v2[0] == 0 && v2[1] == 0)
                return true;
            v2 = v2 / norm(v2);

            if (dot(v2, v0) > dot(v1, v0))
                return false;

            v0 = v1;
        }

        return true;
    };
};

}

#endif /* __CONVEX_IMAGE_HXX__ */
