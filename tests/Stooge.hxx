#ifndef __STOOGE_HXX__
#define __STOOGE_HXX__

#include <iostream>

using namespace std;

class Stooge
{
  public:
    // Factory Method
    static Stooge *make_stooge(int choice);
    virtual void slap_stick() = 0;
};

#endif /* __STOOGE_HXX__ */
