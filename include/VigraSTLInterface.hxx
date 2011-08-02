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

#ifndef __VIGRA_STL_INTERFACE_HXX__
#define __VIGRA_STL_INTERFACE_HXX__

#include <stdarg.h>
#include <vector>
#include <functional>
#include <numeric>
#include "vigra/multi_array.hxx"
#include "objectFeatures.hxx"

namespace bot 
{

/*! A class that provides some interfacing functions for stl data structure 
 *  and vigra data structure
 *
 */
class VigraSTLInterface 
{
public:
	/*! Unfold a vector of matrix by merging the content into a single one vertically
	 *  @param in The input, a std::vector<vigra::MultiArray<2, T > > object
     *  @return The sum of all elements
	 */
    template<class T >
    static vigra::MultiArray<2, T > unfold(const std::vector<vigra::MultiArray<2, T > >& mat_vec)
    {
        vigra::MultiArray<2, T > mat;
        for (int32 ind = 0; ind < mat_vec.size(); ind ++) {
            if (mat.size() == 0)
                mat = mat_vec[ind];
            else
                mat = joinVertically(mat, mat_vec[ind]);
        }

        return mat;
    };

	/*! Compute the sum of all elements from a vector of matrices
	 *  @param in The input, a std::vector<vigra::MultiArray<2, T > > object
     *  @return The sum of all elements
	 */
    template<class T >
    static T accumulate(const std::vector<vigra::MultiArray<2, T > >& in) 
    {
        T sum = 0;
        for (int32 ind = 0; ind < in.size(); ind ++) {
            sum += std::accumulate(
                in[ind].begin(), 
                in[ind].end(), 
                static_cast<T >(0));
        }

    	return sum;
    };

	/*! Convert a std::string to tokens (with given separator) 
	 *  @param str The input string
     *  @param sep The separator
     *  @return A std::vector<str::string > contains the tokens
	 */
    static std::vector<std::string > tokenize(const std::string& str, const char sep=' ') 
    {
    	std::vector<std::string > tokens;
    	std::string tmp = str;
    	size_t pos = tmp.find_first_of(sep);
    	while (pos != std::string::npos) {
    		tokens.push_back(tmp.substr(0, pos));
    		tmp = tmp.substr(pos+1, tmp.length());
    		pos = tmp.find_first_of(sep);
    	}
    	tokens.push_back(tmp.c_str());

    	return tokens;
    };

	/*! Convert a std::string to tokens (with given separator) 
	 *  @param str The input string
     *  @param sep The separator
	 */
    template<class t_elem >
    static std::vector<t_elem > string_to_vector(const std::string& str, const char sep=' ') 
    {
    	std::vector<t_elem > tokens;
    	std::string tmp = str;
    	size_t pos = tmp.find_first_of(sep);
    	while (pos != std::string::npos) {
    		tokens.push_back(static_cast<t_elem >(atof(tmp.substr(0, pos).c_str())));
    		tmp = tmp.substr(pos+1, tmp.length());
    		pos = tmp.find_first_of(sep);
    	}
    	tokens.push_back(static_cast<t_elem >(atof(tmp.c_str())));

    	return tokens;
    };

	/*! Convert a Matrix2D object to a three_set object defined in vigranumpy-extensions
	 *  @param pixels A Matrix2D object
     *  @return A three_set object
	 */
    static three_set pixels2three_set(const Matrix2D& pixels)  
    {
        three_set coords;
        for (int ind = 0; ind < pixels.shape(0); ind ++) {
            if (pixels.shape(1) == 2)
                coords.push_back(three_coordinate(pixels(ind, 0), pixels(ind, 1), 0));
            else
                coords.push_back(three_coordinate(pixels(ind, 0), pixels(ind, 1), pixels(ind, 2)));
        }

        return coords;
    };

	/*! Convert a Matrix2D object to a value_set object defined in vigranumpy-extensions
	 *  @param pixels A Matrix2D object
     *  @return A value_set object
	 */
    static value_set values2value_set(const Matrix2D& values)  
    {
        value_set intens;
        for (int ind = 0; ind < values.shape(0); ind ++) 
            intens.push_back(values(ind, 0));

        return intens;
    };

	/*! Create a vigra matrix by fillding values
	 *  @param count The numbe of values
     *  @param ... Continued values
     *  @return A vigra MultiArray object
	 */
    template<class t_mat_elem >
    static vigra::MultiArray<2, t_mat_elem > fill_matrix(int count, ... )
    {
        if (count <= 0)
            return vigra::MultiArray<2, t_mat_elem >();

        va_list elems;
        va_start(elems, count); 
        typename vigra::MultiArray<2, t_mat_elem >::difference_type shape(1, count);
        vigra::MultiArray<2, t_mat_elem > mat(shape, static_cast<t_mat_elem >(0));
        for (int i = 0; i < count; i ++)
            mat[i] = va_arg(elems, t_mat_elem);
        va_end(elems);

        return mat;
    }

	/*! Fill the values from a std::vector to a vigra::MultiArrayView
	 *  @param vec The std::vector object as the source
     *  @param view The vigra::MultiArrayView object as the target
	 */
    template<class t_mat_elem >
    static void fill_vector_to_view(
        const std::vector<t_mat_elem >& vec, vigra::MultiArrayView<2, t_mat_elem > view) 
    {
        for (int i = 0; i < vec.size(); i ++)
            view[i] = vec[i];
    };

	/*! Fetch the values from a vigra::MultiArrayView and fill them to a std::vector
     *  @param view The vigra::MultiArrayView object as the source
	 *  @return A std::vector object as the target
	 */
    template<class t_mat_elem >
    static std::vector<t_mat_elem > view_to_vector(
        const vigra::MultiArrayView<2, t_mat_elem >& view) 
    {
        std::vector<t_mat_elem > vec;
        for (int i = 0; i < view.size(); i ++)
            vec.push_back(view[i]);

        return vec;
    };

	/*! Fetch the values from a std::vector and fill them to a vigra::MultiArrayView
     *  @param view The std::vector object as the source
	 *  @return A vigra::MultiArray object as the target
	 */
    template<class t_vector_elem, class t_mat_elem >
    static vigra::MultiArray<2, t_mat_elem > vector_to_matrix(
        const std::vector<t_vector_elem >& vec) 
    {
        vigra::MultiArray<2, t_mat_elem > mat;
        VigraSTLInterface::vector_to_matrix<t_vector_elem, t_mat_elem >(vec, mat);

        return mat;
    };

	/*! Fetch the values from a std::vector and fill them to a vigra::MultiArrayView
     *  @param view The std::vector object as the source
	 *  @param mat A vigra::MultiArray object as the target
	 */
    template<class t_vector_elem, class t_mat_elem >
    static void vector_to_matrix(
        const std::vector<t_vector_elem >& vec, 
        vigra::MultiArray<2, t_mat_elem >& mat) 
    {
        int count = vec.size();
        if (count == 0)
            return ;

        int dim = vec[0].size();
        typename vigra::MultiArray<2, t_mat_elem >::difference_type shape(count, dim);
        mat.reshape(shape, static_cast<t_mat_elem >(0));
        for (int i = 0; i < count; i ++) 
            for (int j = 0; j < dim; j ++) 
                mat(i, j) = vec[i][j];
    };

	/*! Fetch the values from a std::vector and fill them to a vigra::MultiArrayView
     *  @param view The std::vector object as the source
	 *  @return A vigra::MultiArray object as the target
	 */
    template<class t_mat_elem >
    static vigra::MultiArray<2, t_mat_elem > vector_to_matrix(
        const std::vector<t_mat_elem >& vec) 
    {
        vigra::MultiArray<2, t_mat_elem > mat;
        VigraSTLInterface::vector_to_matrix<t_mat_elem >(vec, mat);

        return mat;
    };

	/*! Fetch the values from a std::vector and fill them to a vigra::MultiArrayView
     *  @param view The std::vector object as the source
	 *  @param mat A vigra::MultiArray object as the target
	 */
    template<class t_mat_elem >
    static void vector_to_matrix(
        const std::vector<t_mat_elem >& vec, 
        vigra::MultiArray<2, t_mat_elem >& mat) 
    {
        int count = vec.size();
        if (count == 0)
            return ;

        typename vigra::MultiArray<2, t_mat_elem >::difference_type shape(count, 1);
        mat.reshape(shape, static_cast<t_mat_elem >(0));
        for (int i = 0; i < count; i ++) 
            mat(i) = vec[i];
    };

	/*! Fetch the centers from a list of objects
     *  @param objs The input list of objects
	 *  @return A Matrix2D object as the matrix of centers
	 */
    template<class T >
    static Matrix2D get_centers(const std::vector<T >& objs)
    {
        Matrix2D centers;
        for (int ind = 0; ind < objs.size(); ind ++) {
            if (centers.size() == 0) 
                centers.reshape(Matrix2D::difference_type(objs.size(), objs[ind].dim()), 0);

            rowVector(centers, ind) = rowVector(objs[ind].center(), 0);
        }

        return centers;
    };

    //template<class t_obj >
    //static Matrix2D get_centers(const t_obj& objs)
    //{
    //    Matrix2D centers;
    //    for (int ind = 0; ind < objs.size(); ind ++) {
    //        if (centers.size() == 0) 
    //            centers.reshape(Matrix2D::difference_type(objs.size(), objs[ind].dim()), 0);

    //        rowVector(centers, ind) = rowVector(objs[ind].center(), 0);
    //    }

    //    return centers;
    //};
};

}

#endif /* __VIGRA_STL_INTERFACE_HXX__ */
