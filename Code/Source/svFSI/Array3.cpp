
#include "Array3.h"
#include "utils.h"

template<>
double Array3<double>::memory_in_use = 0;

template<>
double Array3<double>::memory_returned = 0;

template<>
int Array3<double>::allocated = 0;

template<>
int Array3<double>::active = 0;

template<>
bool Array3<double>::write_disabled = true;

template<>
void Array3<double>::memory(const std::string& prefix)
{
  utils::print_mem("Array3<double>", prefix, memory_in_use, memory_returned);
}

template<>
void Array3<double>::stats(const std::string& prefix)
{
  utils::print_stats("Array3<double>", prefix, allocated, active);
}

//--------------------------//
//          int             //
//--------------------------//

template<>
double Array3<int>::memory_in_use = 0;

template<>
double Array3<int>::memory_returned = 0;

template<>
int Array3<int>::allocated = 0;

template<>
int Array3<int>::active = 0;

template<>
bool Array3<int>::write_disabled = true;