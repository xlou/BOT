#ifndef __LARRY_HXX__
#define __LARRY_HXX__

#include "Stooge.hxx"

using namespace std;

class Larry: public Stooge
{
  public:
    void slap_stick()
    {
        cout << "Curly: slap head\n";
    }
};

#endif /* __LARRY_HXX__ */
