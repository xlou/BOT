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

#ifndef __EVENT_CONFIGURATION_HXX__
#define __EVENT_CONFIGURATION_HXX__

#include <iostream>
#include <vector>
#include "TypeDefinition.hxx"
#include "SimpleIni/SimpleIni.h"
#include "VigraSTLInterface.hxx"

namespace bot 
{

/*! A class that loads the event definitions from an ini file and manages 
 *  those event definitions.
 *
 */
class EventConfiguration 
{
public:
	/*! Constructor with the input file name specified
	 *  @param filename Path to the event configuration file
	 */
    EventConfiguration(const std::string& filename) 
    {
        load_configuration(filename);
    };

	/*! Load the configuration from a given file 
	 *  @param filename Path to the event configuration file
	 */
    void load_configuration(const std::string& filename)
    {
        CSimpleIniA ini;
        ini.LoadFile(filename.c_str());
        
        CSimpleIniA::TNamesDepend sections;
        ini.GetAllSections(sections);
        sections.sort(CSimpleIniA::Entry::LoadOrder());
        CSimpleIniA::TNamesDepend::const_iterator it;
        for (it = sections.begin(); it != sections.end(); it ++) { 
            std::string section(it->pItem);

            // the global settings ?
            if (section.compare("Global") == 0 || section.compare("global") == 0) {
                const char * pszValue = ini.GetValue(section.c_str(), "k_neighbors", "4");
                k_ = static_cast<int32 >(atof(pszValue));

                pszValue = ini.GetValue(section.c_str(), "distance_threshold", "40");
                d_max_ = static_cast<MatrixElem >(atof(pszValue));

                continue ;
            }

            // default values
            EventPairing pairing;
            MeasureSetting setting;

            // walk through all keys
            CSimpleIniA::TNamesDepend keys;
            ini.GetAllKeys(it->pItem, keys);
            keys.sort(CSimpleIniA::Entry::LoadOrder());
            CSimpleIniA::TNamesDepend::const_iterator it_key;
            std::vector<MatrixElem > w_vec;
            for (it_key = keys.begin(); it_key != keys.end(); it_key ++) { 
                std::string key(it_key->pItem);
                const char * pszValue = ini.GetValue(section.c_str(), key.c_str(), NULL);

                std::transform(key.begin(), key.end(), key.begin(), ::tolower);                
                if (key.find("pairing") != std::string::npos) {
                    pairing = parse_pairing(std::string(pszValue));
                }
                else if (key.find("feature") != std::string::npos) {
                    MatrixElem val;
                    FeatureMeasure measure;
                    parse_joint_feature(std::string(pszValue), measure, val);
                    setting.push_back(measure);
                    w_vec.push_back(val);
                }
                else {
                    std::cerr << "*Warning* unknown key: " << key << std::endl;
                }
            }

            // save this event
            add_event_configuration(
                section, 
                pairing, 
                setting, 
                VigraSTLInterface::vector_to_matrix(w_vec));
        }
    };

	/*! Parse the pairing from an entry from the file
	 *  @param str A std::string object
     *  @return An EventPairing object as a pair of 0 or 1
	 */ 
    EventPairing parse_pairing(const std::string& str)
    {
        std::vector<int32 > tokens = VigraSTLInterface::string_to_vector<int32 >(str);
        if (tokens.size() < 2) {
            std::cerr << "*Warning*: incorrect pairing specficiation: " << str << std::endl;
            return std::make_pair(static_cast<int32 >(1), static_cast<int32 >(1));
        }
        return std::make_pair(tokens[0], tokens[1]);
    };

	/*! Parse the selected features and measures
	 *  @param str A std::string object
     *  @return An MeasureSetting object as a pair of feature and measure
	 */ 
    void parse_joint_feature(
        const std::string& str, 
        FeatureMeasure& measure, 
        MatrixElem& weight)
    {
        std::vector<std::string > tokens = VigraSTLInterface::tokenize(str);
        if (tokens.size() == 1) {
            weight = 1;
            measure = std::make_pair(tokens[0], "MeasureEuclideanDistance");
        }
        else if (tokens.size() == 2) {
            weight = 1;
            measure = std::make_pair(tokens[0], tokens[1]);
        }
        else {
            measure = std::make_pair(tokens[0], tokens[1]);
            weight = static_cast<MatrixElem >(atof(tokens[2].c_str()));
        }
    };

	/*! Add a new event configuration
     *  @param name The name of this event
     *  @param setting The joint feature setting
     *  @param pairing The pairing of this event
     *  @param w The joint feature weights
	 */ 
    void add_event_configuration(
        const std::string& name,
        const EventPairing& pairing,
        const MeasureSetting& setting, 
        const Matrix2D& w) 
    {
        names_.push_back(name);
        settings_.push_back(setting);
        pairings_.push_back(pairing);
        weights_.push_back(w);
    };

	/*! Return the names of all object features
     *  @return A vector of std::string containing the object feature names
	 */ 
    std::vector<std::string > get_feature_names()
    {
        std::vector<std::string > names;
        for (int i = 0; i < settings_.size(); i ++) {
            MeasureSetting setting = settings_[i];
            for (int j = 0; j < setting.size(); j ++) {
                std::string& feature = setting[j].first;
                if (std::find(names.begin(), names.end(), feature) == names.end())
                    names.push_back(feature);
            }
        }

        std::sort(names.begin(), names.end());

        return names;
    };

	/*! Return the names of all measures as the joint features
     *  @return A vector of std::string containing the measures as the joint features
	 */ 
    std::vector<std::string > get_measures()
    {
        std::vector<std::string > measures;
        for (int i = 0; i < settings_.size(); i ++) {
            MeasureSetting setting = settings_[i];
            for (int j = 0; j < setting.size(); j ++) {
                std::string& measure = setting[j].second;
                if (std::find(measures.begin(), measures.end(), measure) == measures.end())
                    measures.push_back(measure);
            }
        }

        std::sort(measures.begin(), measures.end());

        return measures;
    };

	/*! Return the number k of k-nearest-neighbor system
     *  @return An int32 value as the the number k of k-nearest-neighbor system
	 */
    const int32 k() const
    {
        return k_;
    };
 
	/*! Return the number of defined events
     *  @return An int32 value as the the number defined events
	 */
    const int32 count() const
    {
        return names_.size();
    };

	/*! Return the maximally allowed spatial distance in the neighborhood system
     *  @return An MatrixElem value as the maximally allowed spatial distance in the neighborhood system
	 */
    const MatrixElem d_max() const
    {
        return d_max_;
    };
        
	/*! Return a constant reference to the event names
     *  @return A std::vector<std::string > object as the event names
	 */
    const std::vector<std::string >& names() const
    {
        return names_;
    };    
    
	/*! Return a reference to the weights
     *  @return A std::vector<std::string > object as the weights
	 */
    std::vector<Matrix2D >& weights() 
    {
        return weights_;
    };

	/*! Return a constant reference to the feature weights
     *  @return A std::vector<std::string > object as the event names
	 */
    const std::vector<Matrix2D >& weights() const
    {
        return weights_;
    };

	/*! Return a vector of weights all equal to a given scalar
     *  @return A std::vector<std::string > object as the event names
	 */
    std::vector<Matrix2D > weights(const MatrixElem& w) const
    {
        std::vector<Matrix2D > weights = weights_;
        for (int32 ind = 0; ind < weights.size(); ind ++)
            std::fill(weights[ind].begin(), weights[ind].end(), static_cast<MatrixElem >(w));

        return weights;
    };

	/*! Return a constant reference to the event pairings
     *  @return A std::vector<EventPairing > object as the event pairings
	 */
    const std::vector<EventPairing >& pairings() const
    {
        return pairings_;
    };

	/*! Return a constant reference to the settings
     *  @return A std::vector<MeasureSetting > object as the settings
	 */
    const std::vector<MeasureSetting >& settings() const
    {
        return settings_;
    };

	/*! Return a constant reference to the name with given index
     *  @param ind The index
     *  @return A std::string object as the event names
	 */
    const std::string& name(const int32 ind) const
    {
        return names_[ind];
    };
    
	/*! Return a constant reference to the event pairing with given index
     *  @param ind The index
     *  @return A EventPairing object as the event pairings
	 */
    const EventPairing& pairing(const int32 ind) const
    {
        return pairings_[ind];
    };

	/*! Return a constant reference to the setting with given index
     *  @param ind The index
     *  @return A MeasureSetting object as the settings
	 */
    const MeasureSetting& setting(const int32 ind) const
    {
        return settings_[ind];
    };

	/*! Print the event names, features and weights
	 */
    void print() const
    {
        for (int32 ind = 0; ind < names_.size(); ind ++) {
            std::cout << "Event: " << names_[ind] << std::endl;
            for (int32 indF = 0; indF < settings_[ind].size(); indF ++) {
                std::cout << weights_[ind][indF] << ", " << settings_[ind][indF].first << std::endl;
            }
            std::cout << std::endl;
        }
    };
protected:
    std::vector<std::string > names_;
    std::vector<EventPairing > pairings_;
    std::vector<MeasureSetting > settings_;
    std::vector<Matrix2D > weights_;
    int32 k_;
    MatrixElem d_max_;
};

}

#endif /* __EVENT_CONFIGURATION_HXX__ */
