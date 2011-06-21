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

#ifndef __MULTIPLETS_GENERATOR_HXX__
#define __MULTIPLETS_GENERATOR_HXX__

#include "Multiplet.hxx"
#include "Singlet.hxx"
#include "VigraSTLInterface.hxx"
#include "ConvexImage.hxx"
#include "NearestNeighborGenerator.hxx"

namespace bot 
{

/*! A class that generates (hypothetical) multiplets from a list of singlets under some neighborhood constraints
 *
 */
class MultipletsGenerator 
{
public:
    /*! Default constructor
     * 
     */
    MultipletsGenerator() 
    {
        d_max_ = LARGEST_DISTANCE;
        k_ = DEFAULT_NEIGHBORS;
    };

    /*! Constructor with the neighborhood number and distance threshold specified
     *  @param k The neighborhood number
     *  @param d_max The maximally allowed spatial distance
     */
    MultipletsGenerator(const int32 k, const MatrixElem d_max)
    {
        d_max_ = d_max;
        k_ = k;
    };

	/*! Return the maximally allowed spatial distance
	 *  @return A MatrixElem value as the maximally allowed spatial distance
	 */
    MatrixElem d_max() const {
        return d_max_;
    }

	/*! Return the neighborhood number
	 *  @return An int32 value as the neighborhood number
	 */
    int32 k() const {
        return k_;
    }

	/*! Generate multiplets from singlets
	 *  @param img The input raw image
     *  @param seg The input segmentation
     *  @param singlets The list of extracted singlets
     *  @return A list of multiplets
	 */
    Multiplets operator()(const Matrix2D& img, const Matrix2D& seg, const Singlets& singlets) 
    {
        // search for id pairs in a neighborhood (kNN + distance thresholding)
        NearestNeighborGenerator nnGen(k(), d_max(), true);
        
        Matrix2D centers = VigraSTLInterface::get_centers<Singlet >(singlets);
        std::vector<IDPair > id_pairs = nnGen(centers, centers);

        // walk through all id paris and create a party for each of it
        Multiplets multiplets;
        int32 id_party = 0;
        std::vector<IDPair >::const_iterator it;
        for (it = id_pairs.begin(); it != id_pairs.end(); it ++) {
            int32 id1 = it->first;
            int32 id2 = it->second;

            // merge their pixels, values and labels
            Matrix2D pixels_ = joinVertically(singlets[id1].pixels(), singlets[id2].pixels());
            Matrix2D values_ = joinVertically(singlets[id1].values(), singlets[id2].values());
            Matrix2D labels_ = joinVertically(singlets[id1].labels(), singlets[id2].labels());

            // compute the convex hull and the max/min x and y
            Matrix2D convex_hull = ConvexImage::convex_hull(pixels_);
            MatrixElem x_min = convex_hull(argMin(columnVector(convex_hull, 0)), 0);
            MatrixElem y_min = convex_hull(argMin(columnVector(convex_hull, 1)), 1);
            MatrixElem x_max = convex_hull(argMax(columnVector(convex_hull, 0)), 0);
            MatrixElem y_max = convex_hull(argMax(columnVector(convex_hull, 1)), 1);

            // search through the rest pixels in the bounding box
            std::vector<Matrix2D > pixels_vec;
            std::vector<MatrixElem > labels_vec;
            std::vector<MatrixElem > values_vec;
            for (MatrixElem x = x_min; x <= x_max; x ++) {
                for (MatrixElem y = y_min; y <= y_max; y ++) {
                    if (seg(x, y) == singlets[id1].label() || seg(x, y) == singlets[id2].label())
                        continue ;

                    if (!ConvexImage::inside(x, y, convex_hull))
                        continue ;

                    // inside the convex hull but not included in the list
                    pixels_vec.push_back(VigraSTLInterface::fill_matrix<MatrixElem >(2, x, y));
                    values_vec.push_back(img(x, y));
                    //labels_vec.push_back(seg(x, y));
                    labels_vec.push_back(static_cast<MatrixElem >(0));
                }
            }

            // append to the party's coords, values and labels
            pixels_ = joinVertically(pixels_, 
                VigraSTLInterface::vector_to_matrix<Matrix2D, MatrixElem >(pixels_vec));
            values_ = joinVertically(values_, VigraSTLInterface::vector_to_matrix<MatrixElem >(values_vec));
            labels_ = joinVertically(labels_, VigraSTLInterface::vector_to_matrix<MatrixElem >(labels_vec));
            
            // create the party object
            multiplets.push_back(Multiplet(id_party ++, id1, id2, pixels_, values_, labels_));
        }

        return multiplets;
    };

protected:
    MatrixElem d_max_;
    int32 k_;
};

}

#endif /* __MULTIPLETS_GENERATOR_HXX__ */
