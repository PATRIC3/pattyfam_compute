#ifndef _clamp_h
#define _clamp_h

#include <limits>
#include <algorithm>

template <typename T>
T clamp_float(float f)
{
    return static_cast<T>(std::clamp(f,
				     static_cast<float>(std::numeric_limits<T>::min()),
				     static_cast<float>(std::numeric_limits<T>::max())));
				     
}


#endif
