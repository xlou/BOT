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

#ifndef __NEAREST_NEIGHBOR_GENERATOR_HXX__
#define __NEAREST_NEIGHBOR_GENERATOR_HXX__

#include <algorithm>
#include <stdexcept>
#include "TypeDefinition.hxx"
#include "VigraSTLInterface.hxx"
#include "InputOutput.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! The class computes nearest neighbors
 *
 */
class NearestNeighborGenerator 
{
public:
	/*! Default constructor, initialize default values
	 *
	 */
    NearestNeighborGenerator()
    {
        k_ = 3;
        d_max_ = LARGEST_DISTANCE;
        symmetrical_ = false;
    };

	/*! Constructor with k, distance threshold and the symmetry property specified
	 *  @param k The neighborhood number
     *  @param d_max The distance threshold
     *  @param symmetrical Whether the neighborhood is symmetrical
	 */
    NearestNeighborGenerator(
        const int32 k, 
        const MatrixElem d_max, 
        const bool symmetrical)
    {
        k_ = k;
        d_max_ = d_max;
        symmetrical_ = symmetrical;
    };

	/*! Return the neighborhood number
	 *  @return An int32 value as the neighborhood number
	 */
    const int32 k() const
    {
        return k_;
    };

	/*! Return the maximally allowed spatial distance
	 *  @return A MatrixElem value as the maximally allowed spatial distance
	 */
    const MatrixElem d_max() const
    {
        return d_max_;
    };
    
	/*! Return the symmetry property
	 *  @return A bool value that is true when the symmetry property is turned on
	 */
    const bool symmetrical() const
    {
        return symmetrical_;
    };

	/*! Return the kth smallest value in a vector
         *
         *  If k is greater than the vector size, the overall smallest value is returned.
         *  @param vec The input vector
	 *  @return The value of the kth smallest value in a vector
	 */
    template<class t_elem >
    t_elem kth_smallest_value(const std::vector<t_elem >& vec)
    {
	size_t myk;
	size_t vec_size = vec.size();
	if( vec_size < 1 ) {
	    throw std::runtime_error("NearestNeighborGenerator::kth_smallest_value(): vector has to contain at leas one element");
	}
	int n_neighbors = vec.size() - (int) std::count(vec.begin(), vec.end(), LARGEST_DISTANCE);
	int k_min = std::min(k(), n_neighbors);
        std::vector<t_elem > vec_ = vec;
        std::nth_element(vec_.begin(), vec_.begin() + k_min, vec_.end());
        return vec_[k_min];
    };

	/*! Return a list of paired indices that satisfies the neighborhood constraints
     *  @param pts1 A list of points as the source
     *  @param pts2 A list of points as the target
	 *  @return A list of paired indices that satisfies the neighborhood constraints
	 */
    std::vector<IDPair > operator() (const Matrix2D& pts1, const Matrix2D& pts2)
    {
        std::vector<IDPair > id_pairs;
	if (pts1.elementCount() == 0 || pts2.elementCount() == 0)
		return id_pairs;

        // compute the distance matrix: think binomial ¨C (a-b)^2=a^2-2ab+b^2
        Matrix2D::difference_type shape1(pts1.shape(0), 1);
        Matrix2D::difference_type shape2(pts2.shape(0), 1);
        Matrix2D 
            mean1(shape1, 0.0), std1(shape1, 0.0), norm1(shape1, 0.0), 
            mean2(shape2, 0.0), std2(shape2, 0.0), norm2(shape2, 0.0);
        rowStatistics(pts1, mean1, std1, norm1);
        rowStatistics(pts2, mean2, std2, norm2);

        Matrix2D d_mat = sqrt(
            sq(repeatMatrix(norm1, 1, norm2.shape(0))) 
            + sq(repeatMatrix(transpose(norm2), norm1.shape(0), 1))
            - (mmul(pts1, transpose(pts2)) *= 2));

        // only keep the upper right part (if necessary), change the rest to maximum value
        if (symmetrical()) {
            for (int32 x = 0; x < d_mat.shape(0); x ++) 
                for (int32 y = 0; y <= x; y ++) 
                    d_mat(x, y) = LARGEST_DISTANCE;
        }
        
        // determine the k nearest neighbors with d_max distance thresholding
        for (int32 row = 0; row < d_mat.shape(0); row ++) {
		// if having less than k neighbors, no need to compute, add them all
//		if (d_mat.shape(1) - row - 1 <= k())
//			continue;
            Matrix2D::view_type view = rowVector(d_mat, row);
            std::vector<MatrixElem > vec = VigraSTLInterface::view_to_vector<MatrixElem >(view);
            MatrixElem d_max_k = std::min(d_max(), kth_smallest_value(vec));
            std::replace_if(vec.begin(), vec.end(), 
                std::bind2nd(std::greater_equal<MatrixElem >(), d_max_k), LARGEST_DISTANCE);
            VigraSTLInterface::fill_vector_to_view(vec, view);
        }

        // find the index of the neighboring pairs
        for (int32 p = 0; p < d_mat.size(); p ++) {
            if (d_mat[p] >= LARGEST_DISTANCE)
                continue ;

            Matrix2D::difference_type coord = d_mat.scanOrderIndexToCoordinate (p);
            id_pairs.push_back(std::make_pair<int32, int32 >(coord[0], coord[1]));
        };

        return id_pairs;
    };

protected:
    int32 k_;
    MatrixElem d_max_; 
    bool symmetrical_;
};

}

#endif /* __NEAREST_NEIGHBOR_GENERATOR_HXX__ */
