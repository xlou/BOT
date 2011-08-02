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

#ifndef __CONTEXT_HXX__
#define __CONTEXT_HXX__

#include "TypeDefinition.hxx"

namespace bot 
{

/*! A class that contains some global information such as the bounding box and the minimal/maximal intensity.
 *
 */
class Context 
{
public:
    /*! Default constructor
     *
     */
    Context() 
    {
    };

    /*! Constructor that initializes based on a single image
     *  @param image The image as an Matrix2D object
     */
    Context(const Matrix2D& image)
    {
        // compute the bounding box
        bounding_box_.reshape(Matrix2D::difference_type(1, 4));
        bounding_box_[0] = 0;
        bounding_box_[1] = 0;
        bounding_box_[2] = image.shape(0);
        bounding_box_[3] = image.shape(1);

        // compute the max/min intensity
        min_max_.first = *std::min_element(image.begin(), image.end());
        min_max_.second = *std::max_element(image.begin(), image.end());
    };

    /*! Constructor that initializes based on an image sequence
     *  @param images The image sequence as a std::vector<Matrix2D > object
     */
    Context(const std::vector<Matrix2D >& images)
    {
        // compute the bounding box
        bounding_box_.reshape(Matrix2D::difference_type(1, 4));
        bounding_box_[0] = 0;
        bounding_box_[1] = 0;
        bounding_box_[2] = images[0].shape(0);
        bounding_box_[3] = images[0].shape(1);

        // compute the max/min intensity
        min_max_.first = *std::min_element(images[0].begin(), images[0].end());
        min_max_.second = *std::max_element(images[0].begin(), images[0].end());

        // re-compute the max/min intensity
        for (int32 ind = 1; ind < images.size(); ind ++) {
            min_max_.first = std::min(
                min_max_.first,
                *std::min_element(images[ind].begin(), images[ind].end()));
            min_max_.second = std::max(
                min_max_.second,
                *std::max_element(images[ind].begin(), images[ind].end()));
        }
    };

	/*! Return a constant reference to the bounding box object
	 *  @return A Matrix2D object as the bounding box
	 */
    const Matrix2D& bounding_box() const
    {
        return bounding_box_;
    };
    
	/*! Return the globally minimal intensity
	 *  @return An MatrixElem value as the globally minimal intensity
	 */
    const MatrixElem min() const
    {
        return min_max_.first;
    };
    
    /*! Return the globally maximal intensity
	 *  @return An MatrixElem value as the globally maximal intensity
	 */
    const MatrixElem max() const
    {
        return min_max_.second;
    };
protected:
    Matrix2D bounding_box_;
    std::pair<MatrixElem, MatrixElem > min_max_;
};

}

#endif /* __CONTEXT_HXX__ */
