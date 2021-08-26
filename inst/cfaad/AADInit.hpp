/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

// some variables that need to be defined once

#include "AAD.h"
#include <algorithm>

namespace cfaad {

//  AAD.cpp

#if AADMULTIOUT
size_t Node::numAdj = 1;
bool Tape::multi = false;
#endif

Tape globalTape;
thread_local Tape* Number::tape = &globalTape;

} // namespace cfaad
