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

#ifndef __SOLUTION_CODER_HXX__
#define __SOLUTION_CODER_HXX__

#include <iostream>
#include "TypeDefinition.hxx"
#include "Singlet.hxx"
#include "Multiplet.hxx"
#include "InputOutput.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! A class that encodes/decodes the tracking solution
 *
 */
class SolutionCoder
{
public:
    /*! Create a map from the segmentation labels to singlet ids
     *  @param singlets The list of singlets
     *  @return A std::map<int32, int32 > object
     */
    std::map<int32, int32 > map_singlets(const std::vector<Singlet >& singlets)
    {
        std::map<int32, int32 > label_to_id;
        for (int32 ind = 0; ind < singlets.size(); ind ++) 
            label_to_id[singlets[ind].label()] = ind;

        return label_to_id;
    };

    /*! Create a map from the segmentation label pairs to multiplet ids
     *  @param multiplets The list of multiplets
     *  @param singlets The list of singlets
     *  @return A std::map<LabelPair, int32 > object
     */
    std::map<LabelPair, int32 > map_multiplets(
        const std::vector<Multiplet >& multiplets,
        const std::vector<Singlet >& singlets)
    {
        std::map<LabelPair, int32 > label_to_id;
        for (int32 ind = 0; ind < multiplets.size(); ind ++) {
            const std::vector<int32 >& comps = multiplets[ind].components();
            int32 label1 = singlets[comps[0]].label();
            int32 label2 = singlets[comps[1]].label();
            LabelPair key;
            key.first = std::min(label1, label2);
            key.second = std::max(label1, label2);
            label_to_id[key] = ind;
        }

        return label_to_id;
    };

    /*! Make the key for the map from the segmentation labels to objects
     *  @param view A matrix view that contains the segmentation label
     *  @param key An int32 value as the key
     */
    void make_key(const View2D& view, int32& key)
    {
        key = view[0];
    };

    /*! Make the key for the map from the segmentation label pairs to objects
     *  @param view A matrix view that contains the segmentation label pairs
     *  @param key An LabelPair object as the key
     */
    void make_key(const View2D& view, LabelPair& key)
    {
        key = std::make_pair(view[0], view[1]);
    };

    /*! Using the label associations to generate the solution
     *  @param map1 The map (from label(s) to object) for the first frame
     *  @param map2 The map (from label(s) to object) for the second frame
     *  @param hypotheses All the hypotheses
     *  @param key An LabelPair object as the key
     *  @param source The segmentation labels of the source
     *  @param target The associated segmentation labels of the target
     */
    template<class T1, class T2 >
    Matrix2D map_label_to_solution(
        const std::map<T1, int32 >& map1,
        const std::map<T2, int32 >& map2,
        const Hypotheses& hypotheses,
        const Matrix2D& source,
        const Matrix2D& target)
    {
        Matrix2D mat(Shape2D(hypotheses.size(), 1), static_cast<MatrixElem >(0));
        for (int32 ind = 0; ind < source.shape(0); ind ++) {
            int32 id1 = 0;
            if (source[ind, 0] != -1) {
                T1 key1;
                make_key(rowVector(source, ind), key1);
                typename std::map<T1, int32 >::const_iterator it1 = map1.find(key1);
                id1 = it1->second;
            }

            int32 id2 = 0;
            if (target[ind, 0] != -1) {
                T2 key2;
                make_key(rowVector(target, ind), key2);
                typename std::map<T2, int32 >::const_iterator it2 = map2.find(key2);
                id2 = it2->second;
            }

            Hypotheses::const_iterator it = std::find(
                hypotheses.begin(), hypotheses.end(), std::make_pair<int32, int32>(id1, id2));
            if (it != hypotheses.end())
                mat[static_cast<int32 >(it - hypotheses.begin())] = static_cast<MatrixElem >(1);
        }

        return mat;
    };

    /*! Parse the label associations to the solution based on the hypotheses
     *  @param associations The associations that gives assignments between 
     *  segmentation labels
     *  @param events The vector of events
     *  @param singlets1 The singlets from the source
     *  @param singlets2 The singlets from the target
     *  @param multiplets1 The multiplets from the source
     *  @param multiplets2 The multiplets from the target
     *  @param solution The return solution
     */
    bool decode(
        const LabelAssociations& associations, 
        const std::vector<Event >& events,
        const Singlets& singlets1,
        const Singlets& singlets2,
        const Multiplets& multiplets1,
        const Multiplets& multiplets2,
        Solution& solution)
    {
        // create label to singlet/multiplet index maps
        std::map<int32, int32 > map_singlets1 = map_singlets(singlets1);
        std::map<LabelPair, int32 > map_multiplets1 = map_multiplets(multiplets1, singlets1);
        std::map<int32, int32 > map_singlets2 = map_singlets(singlets2);
        std::map<LabelPair, int32 > map_multiplets2 = map_multiplets(multiplets2, singlets2);

        std::vector<std::string > names;
        for (int32 ind = 0; ind < associations.size(); ind ++)
            names.push_back(associations[ind].name);

        for (int32 indE = 0; indE < events.size(); indE ++) {
            const Event& e = events[indE];
            std::vector<std::string >::const_iterator it =
                std::find(names.begin(), names.end(), e.name());
            if (it == names.end()) {
                solution.push_back(Matrix2D(
                    Shape2D(e.hypotheses().size(), 1), 
                    static_cast<MatrixElem >(0)));
                continue;
            }
            int32 indA = static_cast<int32 >(it - names.begin());
            const Matrix2D& source = associations[indA].source;
            const Matrix2D& target = associations[indA].target;
            Matrix2D mat;
            if (e.pairing() == std::make_pair<int32, int32 >(1, 1)) {
                mat = map_label_to_solution(
                    map_singlets1, map_singlets2, 
                    e.hypotheses(), source, target);
            }
            else if (e.pairing() == std::make_pair<int32, int32 >(1, 2)) {
                mat = map_label_to_solution(
                    map_singlets1, map_multiplets2, 
                    e.hypotheses(), source, target);
            }
            else if (e.pairing() == std::make_pair<int32, int32 >(0, 1)) {
                mat = map_label_to_solution(
                    map_singlets1, map_singlets2, 
                    e.hypotheses(), source, target);
            }
            else if (e.pairing() == std::make_pair<int32, int32 >(1, 0)) {
                mat = map_label_to_solution(
                    map_singlets1, map_singlets2, 
                    e.hypotheses(), source, target);
            }
            else if (e.pairing() == std::make_pair<int32, int32 >(2, 1)) {
                mat = map_label_to_solution(
                    map_multiplets1, map_singlets2, 
                    e.hypotheses(), source, target);
            }
            else {
                std::cerr << "*Warning* Skipping unsupported event pairing: " << e.pairing() << std::endl;
                return false ;
            }
            solution.push_back(mat);
        }

        return true;
    };

    /*! Convert a label or label pair to a matrix
     *  @param singlets The list of all singlets
     *  @param object The singlet or multiplet to process
     *  @return A Matrix2D object
     */
    template<class T >
    Matrix2D make_labels(
        const std::vector<Singlet >& singlets,
        const T& object)
    {
        const std::vector<int32 >& components = object.components();
        Matrix2D mat(Shape2D(1, components.size()), static_cast<MatrixElem >(0));
        for (int32 ind = 0; ind < components.size(); ind ++) {
            int32 id = components[ind];
            mat[0, ind] = singlets[id].label();
        }

        return mat;
    };

    /*! Encode the solution into label associations for storage
     *  @param solution The return solution
     *  @param events The vector of events
     *  @param singlets1 The singlets from the source
     *  @param singlets2 The singlets from the target
     *  @param multiplets1 The multiplets from the source
     *  @param multiplets2 The multiplets from the target
     *  @param associations The return associations that gives assignments between 
     *  segmentation labels
     */
    bool encode(
        const Solution& solution, 
        const std::vector<Event >& events,
        const Singlets& singlets1,
        const Singlets& singlets2,
        const Multiplets& multiplets1,
        const Multiplets& multiplets2,
        LabelAssociations& associations)
    {
        for (int32 indE = 0; indE < events.size(); indE ++) {
            const Event& e = events[indE];
            const Matrix2D& x = solution[indE];

            // if the matrix contains only zeros, skip it
            int32 n = std::accumulate(x.begin(), x.end(), static_cast<MatrixElem >(0));
            if (n == 0)
                continue ;
            
            LabelAssociation association;

            // add the name first
            association.name = e.name();

            // walk through the x
            Matrix2D source, target;
            int32 ind = 0;
            for (int32 row = 0; row < x.shape(0); row ++) {
                if (x[row] == 0)
                    continue ;

                int32 id1 = e.hypotheses()[row].first;
                int32 id2 = e.hypotheses()[row].second;
                Matrix2D s_, t_;
                if (e.pairing() == std::make_pair<int32, int32 >(1, 1)) {
                    s_ = make_labels(singlets1, singlets1[id1]);
                    t_ = make_labels(singlets2, singlets2[id2]);
                }
                else if (e.pairing() == std::make_pair<int32, int32 >(1, 2)) {
                    s_ = make_labels(singlets1, singlets1[id1]);
                    t_ = make_labels(singlets2, multiplets2[id2]);
                }
                else if (e.pairing() == std::make_pair<int32, int32 >(0, 1)) {
                    s_ = Matrix2D(Shape2D(1, 1), static_cast<MatrixElem >(-1));
                    t_ = make_labels(singlets2, singlets2[id2]);
                }
                else if (e.pairing() == std::make_pair<int32, int32 >(1, 0)) {
                    s_ = make_labels(singlets1, singlets1[id1]);
                    t_ = Matrix2D(Shape2D(1, 1), static_cast<MatrixElem >(-1));
                }
                else if (e.pairing() == std::make_pair<int32, int32 >(2, 1)) {
                    s_ = make_labels(singlets1, multiplets1[id1]);
                    t_ = make_labels(singlets2, singlets2[id2]);
                }
                else {
                    std::cerr << "*Warning* Skipping unsupported event pairing: " << e.pairing() << std::endl;
                    return false ;
                }

                // initialize matrix source/target, if necessary
                if (source.size() == 0)
                    source.reshape(Shape2D(n, s_.shape(1)), static_cast<MatrixElem >(0));
                
                if (target.size() == 0)
                    target.reshape(Shape2D(n, t_.shape(1)), static_cast<MatrixElem >(0));

                // copy the row vector
                rowVector(source, ind) = s_;
                rowVector(target, ind) = t_;
                ind ++;
            }

            // add to the association
            association.source = source;
            association.target = target;
            
            associations.push_back(association);
        }
        
        //std::cout << association << std::endl;

        return true;
    };
};

}

#endif /* __SOLUTION_CODER_HXX__ */
