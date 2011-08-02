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

#include "InputOutput.hxx"
#include "vigra/matrix.hxx"

using namespace bot;
using namespace vigra::linalg;

std::ostream& operator<< (std::ostream& o, const vigra::MultiArray<2, double >& mat)
{
    for (int x = 0; x < mat.shape(0); x ++) {
        int y;
        for (y = 0; y < mat.shape(1); y ++) {
            o << mat(x, y);
            if (y != mat.shape(1)-1) 
                o << " ";
        }
        if (x != mat.shape(0)-1) 
            o << std::endl;
    }

    return o;
};

std::ostream& operator<< (std::ostream& o, const Object& obj)
{
    o << "id = " << obj.id() << 
        "; components = " << obj.components() << 
        "; center = " << obj.center() << 
        "; label = " << obj.label() << 
        "; count = " << obj.count() << 
        "; dim = " << obj.dim();

    return o;
};

std::ostream& operator<< (std::ostream& o, const Singlet& obj)
{
    o << static_cast<Object >(obj);

    return o;
};

std::ostream& operator<< (std::ostream& o, const Multiplet& obj)
{
    o << static_cast<Object >(obj);

    return o;
};

std::ostream& operator<< (std::ostream& o, const Objects& vec)
{
    return print<Object >(o, vec);
};

std::ostream& operator<< (std::ostream& o, const Singlets& vec)
{
    return print<Singlet >(o, vec);
};

std::ostream& operator<< (std::ostream& o, const Multiplets& vec)
{
    return print<Multiplet >(o, vec);
};

template<class t_elem >
std::ostream& operator<< (std::ostream& o, const std::vector<std::pair<t_elem, t_elem > >& vec)
{
    typename std::vector<std::pair<t_elem, t_elem > >::const_iterator it;
    for (it = vec.begin(); it != vec.end(); it ++) {
        o << it->first << " " << it->second;
        if (it != vec.end() -1 )
            o << std::endl;
    }

    return o;
};

std::ostream& operator << (std::ostream& o, const std::pair<int32, int32 >& p)
{
    o << p.first << " " << p.second;

    return o;
};

std::ostream& operator << (std::ostream& o, const Event& e)
{
    o << "Event: " << e.name() << "; Pairing: " << e.pairing() << std::endl
        << "Hypotheses: " << e.hypotheses().size() << "; Feature: " << e.features().shape() << std::endl;

    return o;
};

std::ostream& operator << (std::ostream& o, const std::vector<std::pair<int32, int32 > >& vec)
{
    return print<std::pair<int32, int32 > >(o, vec);
};

std::ostream& operator << (std::ostream& o, const std::vector<double >& vec)
{
    return print<double >(o, vec);
};

std::ostream& operator << (std::ostream& o, const std::vector<int32 >& vec)
{
    return print<int32 >(o, vec);
};

std::ostream& operator << (std::ostream& o, const std::vector<std::string >& vec)
{
    return print<std::string >(o, vec);
};

std::ostream& operator << (std::ostream& o, const std::vector<vigra::MultiArray<2, double > >& vec)
{
    return print<vigra::MultiArray<2, double > >(o, vec);
};

template<class T >
std::ostream& print(std::ostream& o, const std::vector<T >& vec)
{
    for (int ind = 0; ind < vec.size(); ind ++) {
        o << vec[ind];
        if (ind != vec.size() - 1)
            o << std::endl;
    }
    
    return o;
};

template<class T1, class T2 >
std::ostream& print(std::ostream& o, const std::map<T1, T2 >& m)
{
    typename std::map<T1, T2 >::const_iterator it;
    for (it = m.begin(); it != m.end(); it ++) 
        o << it->first << " - " << it->second << std::endl;
    
    return o;
};

std::ostream& operator<< (std::ostream& o, const EventConfiguration& conf)
{
    const std::vector<std::string >& names = conf.names();
    const std::vector<EventPairing >& pairings = conf.pairings();
    const std::vector<MeasureSetting >& settings = conf.settings();
    const std::vector<Matrix2D >& weights = conf.weights();
    o << "Global: " << std::endl
        << "\tk = " << conf.k() << std::endl
        << "\td_max = " << conf.d_max() << std::endl;

    for (int32 ind = 0; ind < names.size(); ind ++) {
        o << "Event: " << names[ind] << std::endl
            << "\tpairing = " << pairings[ind] << std::endl
            << "\tweights = " << std::endl << weights[ind] << std::endl
            << "\tsetting = " << std::endl << settings[ind] << std::endl;
    }

    return o;
};

std::ostream& operator<< (std::ostream& o, const Context& c)
{
    o << "bounding box: " << c.bounding_box() << std::endl
        << "min - max: " << c.min() << " - " << c.max();

    return o;
};

std::ostream& operator<< (std::ostream& o, const std::map<int32, int32 >& m)
{
    print<int32, int32 >(o, m);
    return o;
};

std::ostream& operator<< (std::ostream& o, const std::map<std::pair<int32, int32 >, int32 >& m)
{
    print<std::pair<int32, int32 >, int32 >(o, m);
    return o;
};

std::ostream& operator<< (std::ostream& o, const LabelAssociations& vec)
{
    print<LabelAssociation >(o, vec);

    return o;
};

std::ostream& operator<< (std::ostream& o, const LabelAssociation& a)
{
    o << a.name << std::endl
        << joinHorizontally(a.source, a.target) << std::endl;

    return o;
};
