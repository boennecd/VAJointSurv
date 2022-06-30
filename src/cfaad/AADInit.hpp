
 /*
 Copyright (c) 2018 Antoine Savine

 This code is from the implementation of the book

 Modern Computational Finance: AAD and Parallel Simulations
 Antoine Savine
 Wiley, 2018

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
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
#ifdef _OPENMP
thread_local Tape* Number::tape = &globalTape;
#else
Tape* Number::tape = &globalTape;
#endif

} // namespace cfaad
