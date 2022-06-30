
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

//  Implementation of the Tape

//  Unchanged for AADET of chapter 15

#include "blocklist.h"
#include "AADNode.h"
#include "AADConfig.h"

namespace cfaad {

constexpr size_t BLOCKSIZE  = 16384;		//	Number of nodes
constexpr size_t ADJSIZE    = 32768;		//	Number of adjoints
constexpr size_t DATASIZE   = 65536;		//	Data in bytes

class Tape
{
#if AADMULTIOUT
	//	Working with multiple results / adjoints?
	static bool							multi;

	//  Storage for adjoints in multi-dimensional case (chapter 14)
    blocklist<double, ADJSIZE>			myAdjointsMulti;
#endif

	//  Storage for derivatives and child adjoint pointers
	blocklist<double, DATASIZE>			myDers;
	blocklist<double*, DATASIZE>		myArgPtrs;

    //  Storage for the nodes
	blocklist<Node, BLOCKSIZE>		    myNodes;

    //  Storage for working memory
    blocklist<double, DATASIZE>         myWKMem;

	//	Padding so tapes in a vector don't interfere
    char                                myPad[64];

#if AADMULTIOUT
    friend auto setNumResultsForAAD(const bool, const size_t);
    friend struct numResultsResetterForAAD;
#endif
	friend class Number;

public:

    //  Build note in place and return a pointer
    //	N : number of childs (arguments)
    template <size_t N>
    Node* recordNode()
    {
        //  Construct the node in place on tape
        Node* node = myNodes.emplace_back(N);

#if AADMULTIOUT
        //  Store and zero the adjoint(s)
        if (multi)
        {
            node->pAdjoints = myAdjointsMulti.emplace_back_multi(Node::numAdj);
            std::fill(node->pAdjoints, node->pAdjoints + Node::numAdj, 0.0);
        }
#endif

      	//	Store the derivatives and child adjoint pointers unless leaf
      	if constexpr(N > 0)
      	{
      		node->pDerivatives = myDers.emplace_back_multi<N>();
      		node->pAdjPtrs = myArgPtrs.emplace_back_multi<N>();

      	}

        return node;
    }

    // create node at run time
    Node* recordNode(size_t const N)
    {
      Node* node = myNodes.emplace_back(N);

#if AADMULTIOUT
      //  Store and zero the adjoint(s)
      if (multi)
      {
        node->pAdjoints = myAdjointsMulti.emplace_back_multi(Node::numAdj);
        std::fill(node->pAdjoints, node->pAdjoints + Node::numAdj, 0.0);
      }
#endif

      if(N)
      {
        node->pDerivatives = myDers.emplace_back_multi(N);
        node->pAdjPtrs = myArgPtrs.emplace_back_multi(N);

      }

      return node;
    }

    // returns the working memory
    double * getWKMem(const size_t N){
        return myWKMem.emplace_back_multi(N);
    }

    //  Reset all adjoints to 0
	void resetAdjoints()
	{
#if AADMULTIOUT
		if (multi)
		{
			myAdjointsMulti.memset(0);
            return;
		}
#endif
        for (Node& node : myNodes)
        {
            node.mAdjoint = 0;
        }
	}

    //  Clear
    void clear()
    {
#if AADMULTIOUT
        myAdjointsMulti.clear();
#endif
		myDers.clear();
		myArgPtrs.clear();
        myNodes.clear();
        myWKMem.clear();
    }

    //  Rewind
    void rewind()
    {

#ifdef _DEBUG

        //  In debug mode, always wipe
        //      makes it easier to identify errors

		clear();

#else
        //  In release mode, rewind and reuse

#if AADMULTIOUT
		if (multi)
		{
			myAdjointsMulti.rewind();
		}
#endif
		myDers.rewind();
		myArgPtrs.rewind();
		myNodes.rewind();
        myWKMem.rewind();

#endif

    }

    //  Set mark
    void mark()
    {
#if AADMULTIOUT
        if (multi)
        {
            myAdjointsMulti.setmark();
        }
#endif
		myDers.setmark();
		myArgPtrs.setmark();
		myNodes.setmark();
        myWKMem.setmark();
    }

    //  Rewind to mark
    void rewindToMark()
    {
#if AADMULTIOUT
        if (multi)
        {
            myAdjointsMulti.rewind_to_mark();
        }
#endif
		myDers.rewind_to_mark();
		myArgPtrs.rewind_to_mark();
		myNodes.rewind_to_mark();
        myWKMem.rewind_to_mark();
    }

    //  Iterators

    using iterator = blocklist<Node, BLOCKSIZE>::iterator;

    auto begin()
    {
        return myNodes.begin();
    }

    auto end()
    {
        return myNodes.end();
    }

    auto markIt()
    {
        return myNodes.mark();
    }

    auto find(Node* node)
    {
        return myNodes.find(node);
    }
};

} // namespace cfaad
