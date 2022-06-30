
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

#pragma once

//  So we can instrument Gaussians like standard math functions
#include "gaussians.h"
#include <memory>
#include "AADConfig.h"
#include "AADNumWrapper.h"
#include "AADVectorFuncs.h"

namespace cfaad {
//  definition of the function to return working memory for the linear algebra
inline double *getLPWKMem(const size_t n){
    return Number::tape->getWKMem(n);
}

//  Routines for multi-dimensional AAD (chapter 14)
//  Set static context for multi-dimensional AAD

//	RAII: reset dimension 1 on destruction
struct numResultsResetterForAAD
{
	~numResultsResetterForAAD()
	{
#if AADMULTIOUT
		Tape::multi = false;
		Node::numAdj = 1;
#endif
	}
};

//  Routine: set dimension and get RAII resetter
inline auto setNumResultsForAAD(const bool multi = false, const size_t numResults = 1)
{
#if AADMULTIOUT
	Tape::multi = multi;
	Node::numAdj = numResults;
	return std::make_unique<numResultsResetterForAAD>();
#endif
}

//  Other utilities

//	Put collection on tape
template <class IT>
inline void putOnTape(IT begin, IT end)
{
    std::for_each(begin, end, [](Number& n) { n.putOnTape(); });
}

//	Convert collection between double and Number or reverse
template<class It1, class It2>
inline void convertCollection(It1 srcBegin, It1 srcEnd, It2 destBegin)
{
    using destType = std::remove_reference_t<decltype(*destBegin)>;
    std::transform(srcBegin, srcEnd, destBegin,
        [](const auto& source) { return destType(source); });
}

} // namespace cfaad
