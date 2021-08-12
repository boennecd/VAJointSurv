#pragma once

//  Use traditional AAD of chapter 10 (false)
//      or expression templated (AADET) of chapter 15 (true)
#ifndef AADET
#define AADET true
#endif

// makes it possible to have multiple derivatives at once
#ifndef AADMULTIOUT
#define AADMULTIOUT false
#endif

// is the library linked with Lapack
#ifndef AADLAPACK
#define AADLAPACK true
#endif

namespace cfaad {
template<class T>
using it_value_type = typename std::iterator_traits<T>::value_type;

template<class T, class U>
using is_it_value_type = std::is_same<it_value_type<T>, U>;
} // namespace cfaad