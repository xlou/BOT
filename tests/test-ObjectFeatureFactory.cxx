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

#include "../include/ObjectFeatureFactory.hxx"
#include "../include/Object.hxx"
#include <vector>

using namespace bot;

int main()
{
    // create object 1 & context 1
    Matrix2D bounding_box1(Matrix2D::difference_type(1, 4));
    bounding_box1[0] = 0;
    bounding_box1[1] = 0;
    bounding_box1[2] = 30;
    bounding_box1[3] = 30;
    Context context1;
    context1.set_bounding_box(bounding_box1);

    Matrix2D pixels1(Matrix2D::difference_type(50, 2));
    Matrix2D values1(Matrix2D::difference_type(50, 1));
    Matrix2D labels1(Matrix2D::difference_type(50, 1));
    int ind = 0;
    for(int x = 0; x < 5; x++) {
        for(int y = 10; y < 20; y++) {
            pixels1(ind, 0) = x;
            pixels1(ind, 1) = y;

            if (x >= 12 && x <= 13 && y >= 23 && y <= 26)
                values1(ind, 0) = 20;
            else
                values1(ind, 0) = 10;

            labels1(ind, 0) = 1;

            ind ++;
        }
    }
    Object obj1(1, pixels1, values1, labels1);

    // create object 2 & context 2
    Matrix2D bounding_box2(Matrix2D::difference_type(1, 6));
    bounding_box2[0] = 0;
    bounding_box2[1] = 0;
    bounding_box2[2] = 0;
    bounding_box2[3] = 100;
    bounding_box2[4] = 300;
    bounding_box2[5] = 500;
    Context context2;
    context2.set_bounding_box(bounding_box2);

    Matrix2D pixels2(Matrix2D::difference_type(75000, 3));
    Matrix2D values2(Matrix2D::difference_type(75000, 1));
    Matrix2D labels2(Matrix2D::difference_type(75000, 1));
    ind = 0;
    for(int x = 0; x < 50; x++) {
        for(int y = 200; y < 225; y++) {
            for(int z = 400; z < 460; z++) {
                pixels2(ind, 0) = x;
                pixels2(ind, 1) = y;
                pixels2(ind, 2) = z;

                if (x >= 100 && x <= 130 && y >= 210 && y <= 225)
                    values2(ind, 0) = 50;
                else
                    values2(ind, 0) = 10;

                labels2(ind, 0) = 5;

                ind ++;
            }
        }
    }
    Object obj2(2, pixels2, values2, labels2);

    // dynamically create feature extractors
    std::string list[] = {
        "ObjectFeatureVolume", "ObjectFeatureSpatialDensityVariation", 
        "ObjectFeaturePosition", "ObjectFeatureBoundingBox", 
        "ObjectFeaturePrincipalComponents", "ObjectFeatureEccentricity", 
        "ObjectFeatureIntensityMean", "ObjectFeatureIntensityDeviation", 
        "ObjectFeatureIntensitySum", "ObjectFeatureIntensityHistogram",         
        "ObjectFeatureBorderDistance", "ObjectFeatureBorderOverlap", 
        "ObjectFeatureNull"
    };
    std::vector<ObjectFeatureFactory* > features;
    for (int i = 0; i < sizeof(list) / sizeof(std::string *); i ++) {
        ObjectFeatureFactory* fea = ObjectFeatureFactory::make(list[i]);
        if (!fea) {
            std::cerr << "Warning: skipping unknown feature class name " << list[i] << std::endl;
        }
        else {
            features.push_back(fea);
        }
    }

    Matrix2D feature_mat;
    for (int i = 0; i < features.size(); i++) {
        std::cout << "**** Extracting " << list[i] << " ****" << std::endl;
        features[i]->extract(feature_mat, obj1, context1);

        std::cout << "Feature " << i << " for object 1: " << std::endl;
        std::cout << "\t\tShape: " << feature_mat.shape() << std::endl;
        std::cout << "\t\tValue: " << feature_mat << std::endl;
    }

    for (int i = 0; i < features.size(); i++)
        delete features[i];

    return 0;
}
