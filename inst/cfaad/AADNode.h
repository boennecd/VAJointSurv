
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