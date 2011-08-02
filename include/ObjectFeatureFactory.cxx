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

#include "ObjectFeatureFactory.hxx"
#include "ObjectFeatureVolume.hxx"
#include "ObjectFeaturePosition.hxx"
#include "ObjectFeatureBoundingBox.hxx"
#include "ObjectFeaturePrincipalComponents.hxx"
#include "ObjectFeatureEccentricity.hxx"
#include "ObjectFeatureIntensityMean.hxx"
#include "ObjectFeatureIntensityDeviation.hxx"
#include "ObjectFeatureIntensitySum.hxx"
#include "ObjectFeatureIntensityHistogram.hxx"
#include "ObjectFeatureSpatialDensityVariation.hxx"
#include "ObjectFeatureBorderDistance.hxx"
#include "ObjectFeatureBorderOverlap.hxx"
#include "ObjectFeatureShapeCompactness.hxx"
#include "ObjectFeatureMassEvenness.hxx"
#include "ObjectFeatureVolumeEvenness.hxx"
#include "ObjectFeatureCenters.hxx"
#include <iostream>
#include <typeinfo>

namespace bot {

ObjectFeatureFactory *ObjectFeatureFactory::make(
    const std::string& classname)
{
    ObjectFeatureFactory *fea = 0;
    if (classname.compare(ObjectFeatureVolume::getClassName()) == 0) {
        fea = new ObjectFeatureVolume();
    }
    else if (classname.compare(ObjectFeaturePosition::getClassName()) == 0) {
        fea = new ObjectFeaturePosition();
    }
    else if (classname.compare(ObjectFeatureBoundingBox::getClassName()) == 0) {
        fea = new ObjectFeatureBoundingBox();
    }
    else if (classname.compare(ObjectFeaturePrincipalComponents::getClassName()) == 0) {
        fea = new ObjectFeaturePrincipalComponents();
    }
    else if (classname.compare(ObjectFeatureEccentricity::getClassName()) == 0) {
        fea = new ObjectFeatureEccentricity();
    }
    else if (classname.compare(ObjectFeatureIntensityMean::getClassName()) == 0) {
        fea = new ObjectFeatureIntensityMean();
    }
    else if (classname.compare(ObjectFeatureIntensityDeviation::getClassName()) == 0) {
        fea = new ObjectFeatureIntensityDeviation();
    }
    else if (classname.compare(ObjectFeatureIntensitySum::getClassName()) == 0) {
        fea = new ObjectFeatureIntensitySum();
    }
    else if (classname.compare(ObjectFeatureIntensityHistogram::getClassName()) == 0) {
        fea = new ObjectFeatureIntensityHistogram();
    }
    else if (classname.compare(ObjectFeatureSpatialDensityVariation::getClassName()) == 0) {
        fea = new ObjectFeatureSpatialDensityVariation();
    }
    else if (classname.compare(ObjectFeatureBorderDistance::getClassName()) == 0) {
        fea = new ObjectFeatureBorderDistance();
    }
    else if (classname.compare(ObjectFeatureBorderOverlap::getClassName()) == 0) {
        fea = new ObjectFeatureBorderOverlap();
    }
    else if (classname.compare(ObjectFeatureShapeCompactness::getClassName()) == 0) {
        fea = new ObjectFeatureShapeCompactness();
    }
    else if (classname.compare(ObjectFeatureMassEvenness::getClassName()) == 0) {
        fea = new ObjectFeatureMassEvenness();
    }
    else if (classname.compare(ObjectFeatureVolumeEvenness::getClassName()) == 0) {
        fea = new ObjectFeatureVolumeEvenness();
    }
    else if (classname.compare(ObjectFeatureCenters::getClassName()) == 0) {
        fea = new ObjectFeatureCenters();
    }

    return fea;
}

ObjectFeatureFactory *ObjectFeatureFactory::make(
    const std::string& classname,
    const Matrix2D& param) 
{
    ObjectFeatureFactory *fea = make(classname);    
    if (fea) {
        fea->set_param(param);
    }

    return fea;
}

}