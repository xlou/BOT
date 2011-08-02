#ifndef __MOE_HXX__
#define __MOE_HXX__

#include "Stooge.hxx"

using namespace std;

class Moe: public Stooge
{
  public:
    void slap_stick()
    {
        cout << "Moe: slap head\n";
    }
};

#endif /* __MOE_HXX__ */
