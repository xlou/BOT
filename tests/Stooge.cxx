#include "Stooge.hxx"
#include "Moe.hxx"
#include "Larry.hxx"
#include "Curly.hxx"
#include <iostream>

using namespace std;

Stooge *Stooge::make_stooge(int choice)
{
  if (choice == 1)
    return new Larry;
  else if (choice == 2)
    return new Moe;
  else
    return new Curly;
}
