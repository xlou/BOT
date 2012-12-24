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

#ifndef __INPUT_OUTPUT_HXX__
#define __INPUT_OUTPUT_HXX__

#include <stdio.h>
#include "Singlet.hxx"
#include "Multiplet.hxx"
#include "Event.hxx"
#include "Context.hxx"
#include "EventConfiguration.hxx"
#include "objectFeatures.hxx"
#include "VigraSTLInterface.hxx"
#include <algorithm>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <dirent.h>
#include "vigra/matrix.hxx"

using namespace bot;

template<class T >
std::ostream& print(std::ostream& o, const std::vector<T >& vec);

template<class T1, class T2 >
std::ostream& print(std::ostream& o, const std::map<T1, T2 >& m);

std::ostream& operator<< (std::ostream& o, const vigra::MultiArray<2, double >& mat);

std::ostream& operator<< (std::ostream& o, const Object& obj);

std::ostream& operator<< (std::ostream& o, const Objects& obj_list);

std::ostream& operator<< (std::ostream& o, const Singlet& obj);

std::ostream& operator<< (std::ostream& o, const Singlets& obj_list);

std::ostream& operator<< (std::ostream& o, const Multiplet& obj);

std::ostream& operator<< (std::ostream& o, const Multiplets& obj_list);

template<class t_elem >
std::ostream& operator<< (std::ostream& o, const std::vector<std::pair<t_elem, t_elem > >& vec);

std::ostream& operator<< (std::ostream& o, const std::map<int32, int32 >& map);

std::ostream& operator<< (std::ostream& o, const std::map<int32, double >& map);

std::ostream& operator << (std::ostream& o, const std::pair<int32, int32 >& p);

std::ostream& operator << (std::ostream& o, const std::vector<std::pair<int32, int32 > >& vec);

std::ostream& operator << (std::ostream& o, const std::vector<double >& vec);

std::ostream& operator << (std::ostream& o, const std::vector<int32 >& vec);

std::ostream& operator << (std::ostream& o, const std::vector<std::string >& vec);

std::ostream& operator << (std::ostream& o, const std::vector<Event >& vec);

std::ostream& operator << (std::ostream& o, const std::vector<vigra::MultiArray<2, double > >& vec);

std::ostream& operator<< (std::ostream& o, const Event& e);

void load_hdf5(const std::string& filename, 
               const std::string& image_path, 
               const std::string& segmentation_path, 
               std::vector<Matrix2D >& images,
               std::vector<Matrix2D >& segmentations);

std::ostream& operator<< (std::ostream& o, const Context& c);

std::ostream& operator<< (std::ostream& o, const EventConfiguration& c);

std::ostream& operator<< (std::ostream& o, const LabelAssociation& a);

std::ostream& operator<< (std::ostream& o, const LabelAssociations& vec);

std::ostream& operator<< (std::ostream& o, const std::map<int32, int32 >& m);

std::ostream& operator<< (std::ostream& o, const std::map<std::pair<int32, int32 >, int32 >& m);


bool stringCompare( const std::string &left, const std::string &right);

std::vector<std::string > getFilesInDir(const std::string& dir, const std::string& filter_ext = ".*");

int countCommonPostfix(const std::string& str1, const std::string& str2);

int findRawFileIndex(const std::vector<std::string >& files, const std::string& file);

// This function parse an annotation string such as "45 -> 43 32"
//bool parseAnnotationLine(const std::string& line, Matrix2D source, Matrix2D target);

#endif /* __INPUT_OUTPUT_HXX__ */
