
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

//  Blocklist data structure for AAD memory management
//  See chapter 10, unchanged with expression templates of chapter 15

#include <array>
#include <list>
#include <iterator>
#include <cstring>
#include <utility>
#include <stdexcept>
#include <memory>

namespace cfaad {

#ifdef DO_CHECKS

/**
 * class used to check that nothing is read after the end of each range. It
 * always allocated with exactly the needed length */
template <class T, size_t block_size>
class blocklist {

  class dum_block {
    size_t n_ele;
    std::unique_ptr<T[]> mem{new T[n_ele]};
    T* begin_v{mem.get()},
     * end_v  {begin_v + n_ele};

  public:
    dum_block(size_t n_ele) : n_ele{n_ele} { }

    T* begin() {
      return begin_v;
    }
    T* begin() const {
      return begin_v;
    }

    T* end() {
      return end_v;
    }
    T* end() const {
      return end_v;
    }

    T& operator[](const size_t i){
      return mem[i];
    }

    const T& operator[](const size_t i) const {
      return mem[i];
    }
  };

  std::list<dum_block>  data;

  using list_iter = decltype(data.begin());
  using block_iter = decltype(data.back().begin());

  list_iter marked_block;

  T* nextblock(const size_t n_ele){
    data.emplace_back(n_ele);
    return data.back().begin();
  }

public:

  blocklist()
  {
    nextblock(1);
  }

  void clear()
  {
    data.clear();
    nextblock(1);
  }

  void rewind()
  {
    clear();
  }

  void memset(unsigned char value = 0)
  {
    for (auto& arr : data)
    {
      std::memset(&arr[0], value,
                  std::distance(arr.begin(), arr.end()) * sizeof(T));
    }
  }

  template<typename ...Args>
  T* emplace_back(Args&& ...args)
  {
    block_iter next_space = nextblock(1);
    T* emplaced = new (&*next_space)T(std::forward<Args>(args)...);
    return emplaced;
  }

  T* emplace_back()
  {
    return nextblock(1);
  }

  template <size_t n>
  T* emplace_back_multi()
  {
    return nextblock(n);
  }

  T* emplace_back_multi(const size_t n)
  {
    return nextblock(n);
  }

  void setmark()
  {
    marked_block = std::prev(data.end());
  }

  //  Rewind to mark
  void rewind_to_mark()
  {
    list_iter marked_block_p1{marked_block};
    ++marked_block_p1;
    if(marked_block != data.end())
      data.erase(marked_block_p1, data.end());
  }

  class iterator
  {
    //  List and block
    list_iter        cur_block;		//  current block
    block_iter       cur_space;		//  current space
    block_iter       first_space;	//  first space in block
    block_iter       last_space;	//  last (+1) space in block

  public:

    //  iterator traits
    using difference_type = ptrdiff_t;
    using reference = T&;
    using pointer = T*;
    using value_type = T;
    using iterator_category = std::bidirectional_iterator_tag;

    iterator() {}
    iterator(list_iter cb, block_iter cs, block_iter fs, block_iter ls) :
      cur_block(cb), cur_space(cs), first_space(fs), last_space(ls) {}

    iterator& operator++()
    {
      if (cur_space == last_space || ++cur_space == last_space)
      {
        ++cur_block;
        first_space = cur_block->begin();
        last_space = cur_block->end();
        cur_space = first_space;
      }

      return *this;
    }

    iterator& operator--()
    {
      if (cur_space == first_space)
      {
        --cur_block;
        first_space = cur_block->begin();
        last_space = cur_block->end();
        cur_space = last_space;
      }

      --cur_space;

      return *this;
    }

    T& operator*()
    {
      return *cur_space;
    }
    const T& operator*() const
    {
      return *cur_space;
    }
    T* operator->()
    {
      return &*cur_space;
    }
    const T* operator->() const
    {
      return &*cur_space;
    }

    //	Check equality
    bool operator ==(const iterator& rhs) const
    {
      return (cur_block == rhs.cur_block && cur_space == rhs.cur_space);
    }
    bool operator !=(const iterator& rhs) const
    {
      return (cur_block != rhs.cur_block|| cur_space != rhs.cur_space);
    }
  };

  iterator begin()
  {
    return iterator(data.begin(), data.begin()->begin(),
                    data.begin()->begin(), data.begin()->end());
  }

  iterator end()
  {
    list_iter cur_block{std::prev(data.end())};
    return iterator(cur_block, cur_block->end(),
                    cur_block->begin(), cur_block->end());
  }

  iterator mark()
  {
    return iterator(marked_block, marked_block->end(),
                    marked_block->begin(), marked_block->end());
  }

  iterator find(const T* const element)
  {
    iterator it = end();
    iterator b = begin();

    while (it != b)
    {
      --it;
      if (&*it == element) return it;
    }

    if (&*it == element) return it;

    return end();
  }
};

#else // #ifdef DO_CHECKS

template <class T, size_t block_size>
class blocklist
{
    //  Container = list of blocks
    std::list<std::array<T, block_size> >  data;

    using list_iter = decltype(data.begin());
    using block_iter = decltype(data.back().begin());

    //  Current block
    list_iter           cur_block;

	//  Last block
	list_iter			last_block;

	//  Next free space in current block
    block_iter          next_space;

    //  Last free space (+1) in current block
    block_iter          last_space;

    //  Mark
    list_iter           marked_block;
    block_iter          marked_space;

    //  Create new array
    void newblock()
    {
        data.emplace_back();
        cur_block = last_block = std::prev(data.end());
        next_space = cur_block->begin();
        last_space = cur_block->end();
    }

    //  Move on to next array
    void nextblock()
    {
        //  This is the last array: create new
        if (cur_block == last_block)
        {
            newblock();
        }
        else
        {
            ++cur_block;
            next_space = cur_block->begin();
            last_space = cur_block->end();
        }
    }

public:

    //  Create first block on construction
    blocklist()
    {
        newblock();
    }

    //  Factory reset
    void clear()
    {
        data.clear();
        newblock();
    }

    //  Rewind but keep all blocks
    void rewind()
    {
        cur_block = data.begin();
        next_space = cur_block->begin();
        last_space = cur_block->end();
    }

	//	Memset
	void memset(unsigned char value = 0)
	{
		for (auto& arr : data)
		{
			std::memset(&arr[0], value, block_size * sizeof(T));
		}
	}

    //  Construct object of type T in place
    //      in the next free space and return a pointer on it
    //  Implements perfect forwarding of constructor arguments
    template<typename ...Args>
    T* emplace_back(Args&& ...args)
    {
        //  No more space in current array
        if (next_space == last_space)
        {
            nextblock();
        }
        //  Placement new, construct in memory pointed by next
        T* emplaced = new (&*next_space)    //  memory pointed by next as T*
            T(std::forward<Args>(args)...);      //  perfect forwarding of ctor arguments

        //  Advance next
        ++next_space;

        //  Return
        return emplaced;
    }

    //  Overload for default constructed
    T* emplace_back()
    {
        //  No more space in current array
        if (next_space == last_space)
        {
            nextblock();
        }

        //  Current space
        auto old_next = next_space;

        //  Advance next
        ++next_space;

        //  Return
        return &*old_next;
    }

	//  Stores n default constructed elements
    //      and returns a pointer on the first

	//	Version 1: n known at compile time
	template <size_t n>
	T* emplace_back_multi()
	{
        static_assert(n <= block_size,
                      "requested number of elements is greater than the block size");

		//  No more space in current array
		if (std::distance(next_space, last_space) <  static_cast
                <typename std::iterator_traits<block_iter>::difference_type>(n))
		{
			nextblock();
		}

		//  Current space
		auto old_next = next_space;

		//  Advance next
		next_space += n;

		//  Return
		return &*old_next;
	}

	//	Version 2: n unknown at compile time
	T* emplace_back_multi(const size_t n)
	{
		//  No more space in current array
		if (std::distance(next_space, last_space) < static_cast
                <typename std::iterator_traits<block_iter>::difference_type>(n))
		{
            if(n > block_size)
                throw std::runtime_error
                    ("requested number of elements is greater than the block size");
			nextblock();
		}

		//  Current space
		auto old_next = next_space;

		//  Advance next
		next_space += n;

		//  Return
		return &*old_next;
	}

	//	Marks

    //  Set mark
    void setmark()
    {
        if (next_space == last_space)
        {
            nextblock();
        }

        marked_block = cur_block;
        marked_space = next_space;
    }

    //  Rewind to mark
    void rewind_to_mark()
    {
        cur_block = marked_block;
        next_space = marked_space;
		last_space = cur_block->end();
    }

    //  Iterator

    class iterator
    {
        //  List and block
        list_iter        cur_block;		//  current block
        block_iter       cur_space;		//  current space
        block_iter       first_space;	//  first space in block
        block_iter       last_space;	//  last (+1) space in block

    public:

        //  iterator traits
        using difference_type = ptrdiff_t;
        using reference = T&;
        using pointer = T*;
        using value_type = T;
        using iterator_category = std::bidirectional_iterator_tag;

        //	Default constructor
        iterator() {}

        //	Constructor
        iterator(list_iter cb, block_iter cs, block_iter fs, block_iter ls) :
            cur_block(cb), cur_space(cs), first_space(fs), last_space(ls) {}

        //	Pre-increment (we do not provide post)
        iterator& operator++()
        {
            ++cur_space;
            if (cur_space == last_space)
            {
                ++cur_block;
                first_space = cur_block->begin();
                last_space = cur_block->end();
				cur_space = first_space;
            }

            return *this;
        }

        //	Pre-decrement
        iterator& operator--()
        {
            if (cur_space == first_space)
            {
                --cur_block;
                first_space = cur_block->begin();
                last_space = cur_block->end();
				cur_space = last_space;
            }

            --cur_space;

            return *this;
        }

        //	Access to contained elements
        T& operator*()
        {
            return *cur_space;
        }
        const T& operator*() const
        {
            return *cur_space;
        }
        T* operator->()
        {
            return &*cur_space;
        }
        const T* operator->() const
        {
            return &*cur_space;
        }

        //	Check equality
        bool operator ==(const iterator& rhs) const
        {
            return (cur_block == rhs.cur_block && cur_space == rhs.cur_space);
        }
        bool operator !=(const iterator& rhs) const
        {
            return (cur_block != rhs.cur_block|| cur_space != rhs.cur_space);
        }
    };

    //  Access to iterators

    iterator begin()
    {
        return iterator(data.begin(), data.begin()->begin(),
            data.begin()->begin(), data.begin()->end());
    }

    iterator end()
    {
        return iterator(cur_block, next_space,
            cur_block->begin(), cur_block->end());
    }

    //  Iterator on mark
    iterator mark()
    {
        return iterator(marked_block, marked_space,
            marked_block->begin(), marked_block->end());
    }

    //  Find element, by pointer, searching sequentially from the end
    iterator find(const T* const element)
    {
        //	Search from the end
        iterator it = end();
        iterator b = begin();

        while (it != b)
        {
            --it;
            if (&*it == element) return it;
        }

        if (&*it == element) return it;

        return end();
    }
};

#endif // #ifdef DO_CHECKS

} // namespace cfaad
