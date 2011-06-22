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

#ifndef __OBJECT_HXX__
#define __OBJECT_HXX__

#include <vector>
#include <vigra/matrix.hxx>
#include "TypeDefinition.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! This class is the base class for Singlet and Multiplet. It contains 
 *  all the basic informaton of an object (e.g. its pixel coordinates, 
 *  intensities, components).
 *
 */
class Object 
{
public:
    /* Default constructor
     *
     */
    Object() 
    {
        id_ = -1;
        components_.push_back(id_);
    };

	/*! Constructor with int32, pixels, labels and values specified
     *  @param id The id of this multiplet
     *  @param pixels The pixels (coordinates) that this multiplet contains
     *  @param values The values (intensities) that this multiplet contains
     *  @param labels The segmentation labels that this multiplet contains
	 */
    Object(const int32 id, const Matrix2D& pixels, const Matrix2D& values, const Matrix2D& labels)
        : id_(id), pixels_(pixels), values_(values), labels_(labels)
    {
        // compute its center
        center_.reshape(Matrix2D::difference_type(1, dim()), static_cast<MatrixElem >(0));
        columnStatistics(pixels_, center_);
    };

	/*! Add a new feature with its name
     *  @param name The name of this feature
     *  @param feature The feature, a Matrix2D object
	 */
    void add_feature(const std::string& name, const Matrix2D& feature) 
    {
        names_.push_back(name);
        features_.push_back(feature);
    };

	/*! Get a feature by its name
     *  @param name The name of this feature
     *  @return The feature, a Matrix2D object
	 */
    const Matrix2D& get_feature(const std::string& name) const
    {
        std::vector<std::string >::const_iterator it = 
            std::lower_bound(names_.begin(), names_.end(), name);
        
        return get_feature(static_cast<int32 >(it - names_.begin()));
    };

	/*! Get a feature by its index
     *  @param ind The ind of this feature
     *  @return The feature, a Matrix2D object
	 */
    const Matrix2D& get_feature(const int32 ind) const
    {
        if (ind < 0)
            return features_.front();
        else if (ind >= features_.size())
            return features_.back();
        
        return features_[ind];
    };

    /*! Return a const reference to the pixel values (intensities)
        @return An Matrix2D type value representing the pixel values
     */
    const Matrix2D& values() const
    {
        return values_;
    };

    /*! Return a const reference to the pixels (coordinates)
        @return An Matrix2D type value representing the pixels
     */
    const Matrix2D& pixels() const
    {
        return pixels_;
    };

    /*! Return a const reference to the labels (segmentation id)
        @return An Matrix2D type value representing the labels
     */
    const Matrix2D& labels() const
    {
        return labels_;
    };

    /*! Return a const reference to the pixels (coordinates)
        @return An Matrix2D type value representing the pixels
     */
    const MatrixElem label() const
    {
        if (labels_.size() == 0)
            return 0;

        return labels_[0];
    };

    /*! Return a const reference to the cener (coordinates)
        @return An Matrix2D type value representing the pixels
     */
    const Matrix2D& center() const
    {
        return center_;
    };

    /*! Return the id of the object
        @return An int32 type value representing the id
     */
    const int32 id() const 
    {
        return id_;
    };

    /*! Return the dimension of the image where the point cloud is extracted
        @return An int32 type value representing the number of dimensions
     */
    const int32 dim() const 
    {
        return pixels_.shape(1);
    };

    /*! Return the number of the points
        @return An int32 type value representing the number of dimensions
     */
    const int32 count() const 
    {
        return pixels_.shape(0);
    };
    
    /*! Return a const reference to the id vector of the components
        @return The components
     */
    const std::vector<int32 >& components() const 
    {
        return components_;
    };

    /*! Return if the object is a singlet
        @return An bool value that is true if this object is a singlet
     */
    const bool is_singlet() const 
    {
        return components_.size() <= 1;
    };

    /* Returne a reference to the feature vector
     * @return A reference to the feature vector
     */
    const std::vector<Matrix2D >& features() const
    {
        return features_;
    };

    /* Returne a reference to the name vector
     * @return A reference to the name vector
     */
    const std::vector<std::string >& names() const
    {
        return names_;
    };
protected:
    /* id of this object */
    int32 id_;

    /* id of the objects this party contains */
    std::vector<int32 > components_;

    /* Coordinate of the center of this object */
    Matrix2D center_;

    /* 2D matrix of the coordinates of the pixels */
    Matrix2D pixels_;

    /* 2D matrix of the intensity values of the pixels */
    Matrix2D values_;

    /* 2D matrix of the labels of the pixels */
    Matrix2D labels_;

    /* Object feature names */
    std::vector<std::string > names_;

    /* Object features */
    std::vector<Matrix2D > features_;
};

typedef std::vector<Object > Objects;
typedef std::vector<Objects > ObjectsSequence;

}

#endif /* __OBJECT_HXX__ */
