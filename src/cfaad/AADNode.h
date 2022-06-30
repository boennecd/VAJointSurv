
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

//  AAD implementation of chapter 10
//  (With multi-dimensional additions of chapter 14)

//  Implementation of Node = record on tape

//  Unchanged for AADET of chapter 15

#include <exception>
#include <algorithm>

namespace cfaad {

class Node
{
	friend class Tape;
	friend class Number;
#if AADMULTIOUT
	friend auto setNumResultsForAAD(const bool, const size_t);
	friend struct numResultsResetterForAAD;
#endif

    //  The adjoint(s)
	//	in single case, self held (chapter 10)
	double			mAdjoint = 0;
#if AADMULTIOUT
	//	in multi case, held separately and accessed by pointer (chapter 14)
    double*         pAdjoints;
#endif

	//  Data lives in separate memory

    //  the n derivatives to arguments,
    double*         pDerivatives;

    //  the n pointers to the adjoints of arguments
    double**        pAdjPtrs;

#if AADMULTIOUT
    //  Number of adjoints (results) to propagate, usually 1
    //  See chapter 14
    static size_t   numAdj;
#endif

    //  Number of childs (arguments)
    const size_t    n;

public:

    Node(const size_t N = 0) : n(N) {}

    //  Access to adjoint(s)
	//	single
    double& adjoint()
    {
		return mAdjoint;
	}
#if AADMULTIOUT
	//	multi
	double& adjoint(const size_t n) { return pAdjoints[n]; }
#endif

    //  Back-propagate adjoints to arguments adjoints

    //  Single case, chapter 10
    void propagateOne()
    {
		//  Nothing to propagate
		if (!n || !mAdjoint) return;

		for (size_t i = 0; i < n; ++i)
        {
			*(pAdjPtrs[i]) += pDerivatives[i] * mAdjoint;
        }
    }

#if AADMULTIOUT
    //  Multi case, chapter 14
    void propagateAll()
    {
        //  No adjoint to propagate
        if (!n || std::all_of(pAdjoints, pAdjoints + numAdj,
            [](const double& x) { return !x; }))
            return;

        for (size_t i = 0; i < n; ++i)
        {
            double *adjPtrs = pAdjPtrs[i], ders = pDerivatives[i];

            //  Vectorized!
            for (size_t j = 0; j < numAdj; ++j)
            {
                adjPtrs[j] += ders * pAdjoints[j];
            }
        }
    }
#endif
};

} // namespace cfaad
