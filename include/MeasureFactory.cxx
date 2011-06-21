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

#include "MeasureFactory.hxx"
#include "MeasureL1L2Distance.hxx"
#include "MeasureEarthMoversDistance.hxx"
#include "MeasureAnglePattern.hxx"
#include "MeasureAppointment.hxx"
#include <iostream>
#include <typeinfo>

namespace bot {

MeasureFactory *MeasureFactory::make(
    const std::string& classname)
{
    MeasureFactory *fea = 0;
    if (classname.compare(MeasureEarthMoversDistance::getClassName()) == 0) {
        fea = new MeasureEarthMoversDistance();
    }
    else if (classname.compare(MeasureAnglePattern::getClassName()) == 0) {
        fea = new MeasureAnglePattern();
    }
    else if (classname.compare(MeasureAppointmentLeft::getClassName()) == 0) {
        fea = new MeasureAppointmentLeft();
    }
    else if (classname.compare(MeasureAppointmentRight::getClassName()) == 0) {
        fea = new MeasureAppointmentRight();
    }
    else if (classname.compare(MeasureEuclideanDistance::getClassName()) == 0) {
        fea = new MeasureEuclideanDistance();
    }
    else if (classname.compare(MeasureNormalizedEuclideanDistance::getClassName()) == 0) {
        fea = new MeasureNormalizedEuclideanDistance();
    }
    else if (classname.compare(MeasureManhattanDistance::getClassName()) == 0) {
        fea = new MeasureManhattanDistance();
    }
    else if (classname.compare(MeasureNormalizedManhattanDistance::getClassName()) == 0) {
        fea = new MeasureNormalizedManhattanDistance();
    }

    return fea;
}

MeasureFactory *MeasureFactory::make(
    const std::string& classname,
    const Matrix2D& param) 
{
    MeasureFactory *fea = make(classname);    
    if (fea) {
        fea->set_param(param);
    }

    return fea;
}

}