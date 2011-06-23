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

#ifndef __SINGLETS_GENERATOR_HXX__
#define __SINGLETS_GENERATOR_HXX__

#include "vigra/matrix.hxx"
#include "Singlet.hxx"
#include "TypeDefinition.hxx"
#include "VigraSTLInterface.hxx"

namespace bot 
{

/*! A class that generates singlets
 *
 */
class SingletsGenerator 
{
public:
    /*! Default constructor
     * 
     */
    SingletsGenerator() 
    {
    };

	/*! Generate singlets from the raw image and the segmentation
	 *  @param img The input raw image
     *  @param seg The input segmentation
     *  @return A list of singlets
	 */
    Singlets operator()(const Matrix2D& img, const Matrix2D& seg) 
    {
        // create several maps: mapping a label in the segmentation to an id, 
        // an coordinate vector, an intensity vector and a label vector        
        std::map<MatrixElem, int32 > seg_to_id;
        std::map<MatrixElem, std::vector<Matrix2D::difference_type > > seg_to_pixel;
        std::map<MatrixElem, std::vector<MatrixElem > > seg_to_image;
        std::map<MatrixElem, std::vector<MatrixElem > > seg_to_label;

        // walk through the entire image pixel by pixel
        int32 base_id = 0;
        for (int p = 0; p < img.elementCount(); p ++) {
            MatrixElem label = seg[p];

            // is background ?
            if (label == 0)
                continue;

            // is there an id assigned to this label already?
            if (seg_to_id.find(label) == seg_to_id.end()) {     // no, then add a new one
                seg_to_id[label] = base_id ++;
                seg_to_pixel[label] = std::vector<Matrix2D::difference_type >();
                seg_to_image[label] = std::vector<MatrixElem >();
                seg_to_label[label] = std::vector<MatrixElem >();
            }

            // append the coord, intensity and label
            seg_to_pixel[label].push_back(seg.scanOrderIndexToCoordinate(p));
            seg_to_image[label].push_back(img[p]);
            seg_to_label[label].push_back(label);

        }

        // create the object list
        Singlets singlets(seg_to_id.size());
        std::map<MatrixElem, int32 >::iterator it;
        for (it = seg_to_id.begin(); it != seg_to_id.end(); it ++) {
            MatrixElem label = it->first;
            int32 id = it->second;

            singlets[id] = Singlet(id, 
                VigraSTLInterface::vector_to_matrix<Matrix2D::difference_type, MatrixElem>(seg_to_pixel[label]),
                VigraSTLInterface::vector_to_matrix<MatrixElem>(seg_to_image[label]),
                VigraSTLInterface::vector_to_matrix<MatrixElem>(seg_to_label[label]));
        }

        return singlets;
    };
};

}

#endif /* __SINGLETS_GENERATOR_HXX__ */
