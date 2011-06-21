#ifndef __CURLY_HXX__
#define __CURLY_HXX__

#include "Stooge.hxx"

using namespace std;

class Curly: public Stooge
{
  public:
    void slap_stick()
    {
        cout << "Curly: slap head\n";
    }
};

#endif /* __CURLY_HXX__ */
