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

#ifndef __EVENT_HXX__
#define __EVENT_HXX__

#include "TypeDefinition.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that contains the information of a defined event, such as name, 
 *  selected features and pairing (e.g. one to many, many to one).
 *
 */
class Event {
public:
    /*! Default constructor
     *
     */
    Event() {};

    /*! Constructor with name, hypotheses and feature matrix specified
     *  @param name The event name
     *  @param pairing The pairing of this event
     *  @param hypotheses The hypotheses of associations
     *  @param features The features for those hypotheses
     */
    Event(
        const std::string& name, 
        const EventPairing& pairing,
        const Hypotheses& hypotheses,
        const Matrix2D& features)
    {
        name_ = name;
        pairing_ = pairing;
        hypotheses_ = hypotheses;
        features_ = features;
    };

	/*! Return a constant reference to the name
	 *  @return A std::string object as the name
	 */
    const std::string& name() const
    {
        return name_;
    };

	/*! Return a constant reference to the pairing
	 *  @return A std::string object as the pairing
	 */
    const EventPairing& pairing() const
    {
        return pairing_;
    };

	/*! Return a constant reference to the hypotheses
	 *  @return A Hypotheses object as the v
	 */
    const Hypotheses& hypotheses() const
    {
        return hypotheses_;
    };

	/*! Return the size of this event, i.e. the number of hypotheses
	 *  @return A int32 values as the size
	 */
    const int32 size() const
    {
        return hypotheses_.size();
    };

	/*! Return a constant reference to the features
	 *  @return An Matrix2D object as the feature matrix
	 */
    const Matrix2D& features() const
    {
        return features_;
    };

	/*! Return a constant reference to the hypotheses
	 *  @return A Hypotheses object as the v
	 */
    Matrix2D& features()
    {
        return features_;
    };
protected:
    std::string name_;
    EventPairing pairing_;
    Hypotheses hypotheses_;
    Matrix2D features_;
};

}

#endif /* __EVENT_HXX__ */
