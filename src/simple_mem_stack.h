#ifndef SIMPLE_MEM_STACK_H
#define SIMPLE_MEM_STACK_H

#include <vector>
#include <iterator>
#include <stack>
#include <list>
#include <cstddef>

#ifdef DEBUG_SIMPLE_MEM_STACK
#include <memory>
#include <algorithm>
#endif

namespace ghqCpp {

/**
 * stack like object used to avoid repeated allocation. In principle,
 * everything goes well if set_mark_raii() is called after all requests in a
 * function call and before any further function calls and subsequently the
 * returned object is destroyed at the end of the scope (as typical with RAII).
 *
 * If you define DEBUG_SIMPLE_MEM_STACK then just the right memory for each
 * request is allocated. Thus, you will easily find errors due to too little
 * memory being requested.
 */
template<class T>
class simple_mem_stack {
#ifdef DEBUG_SIMPLE_MEM_STACK
  class block_container {
    size_t const n_ele;
    std::unique_ptr<T[]> mem{new T[n_ele]};

  public:
    using iterator = T*;
    using const_iterator = const T*;

    block_container(size_t const n_ele): n_ele{n_ele} { }
    block_container(const block_container &o):
      n_ele{o.n_ele},
      mem{new T[n_ele]} {
        std::copy(o.begin(), o.end(), begin());
      }

    iterator begin() { return mem.get(); }
    iterator end() { return mem.get() + n_ele; }
    const_iterator begin() const  { return mem.get(); }
    const_iterator end() const { return mem.get() + n_ele; }
    size_t size() const { return n_ele; }
  };

#else // #ifdef DEBUG_SIMPLE_MEM_STACK
  // TODO: maybe replace with a simpler container?
  using block_container = std::vector<T>;

#endif

  using outer_container = std::list<block_container>;
  /// holds the allocated memory in blocks
  outer_container memory;

  /**
   * a simple iterator that implements the members we need. It also stores an
   * iterator to the block which is needed when go back to a mark.
   */
  class iterator {
    using block_it = typename outer_container::iterator;
    using block_cont_it = typename block_container::iterator;
    block_cont_it cur_ptr;

  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = T;
    using pointer           = T*;
    using reference         = T&;

    block_it cont;

    iterator() = default;
    iterator(block_it cont, block_cont_it cur_ptr):
      cur_ptr{cur_ptr}, cont{cont} { }
    iterator(block_it cont):
      iterator(cont, cont->begin()) { }

    reference operator*() const { return *cur_ptr; }
    pointer operator->() { return cur_ptr; }

    iterator& operator++() {
      cur_ptr++;
      return *this;
    }
    iterator operator++(int) {
      iterator res = *this;
      this->operator++();
      return res;
    }

    iterator& operator+=(const difference_type n){
      cur_ptr += n;
      return *this;
    }
    friend iterator operator+
      (const iterator &it, const difference_type n){
      iterator out{it};
      return out += n;
    }
    friend iterator operator+
      (const difference_type n, const iterator &it){
      return it + n;
    }

    friend bool operator==(const iterator& a, const block_cont_it& b) {
      return a.cur_ptr == b;
    };
    friend bool operator==(const block_cont_it& b, const iterator& a) {
      return a == b;
    };
    friend bool operator==(const iterator& a, const iterator& b) {
      return a == b.cur_ptr;
    };

    friend bool operator!=(const iterator& a, const block_cont_it& b) {
      return a.cur_ptr != b;
    };
    friend bool operator!=(const block_cont_it& b, const iterator& a) {
      return a != b;
    };
    friend bool operator!=(const iterator& a, const iterator& b) {
      return a != b.cur_ptr;
    };

    friend bool operator>=(const iterator& a, const block_cont_it& b) {
      return a.cur_ptr >= b;
    };
    friend bool operator>=(const block_cont_it& b, const iterator& a) {
      return b >= a.cur_ptr;
    };
    friend bool operator>=(const iterator& a, const iterator& b) {
      return a >= b.cur_ptr;
    };
  };

  /// markers to jump back to
  std::stack<iterator> marks;

  /// the current head (next free object). May be an end pointer
  iterator cur_head;

  /// the minimum size of the blocks
#ifdef DEBUG_SIMPLE_MEM_STACK
  static constexpr size_t min_block_size{0};
#else
  static constexpr size_t min_block_size{32768};
#endif

  /**
   * creates a new block of given minimum size. The new block size will be
   * at least as great as n_ele
   */
  void new_block(size_t const n_ele){
    auto it = cur_head.cont;
    while(++it != memory.end() && it->size() < n_ele) { }

    if(it == memory.end()){
      // did not find a block of the appropriate size. Create a new block
#ifdef DEBUG_SIMPLE_MEM_STACK
      // allocate just enough memory
      memory.emplace_back(n_ele);
#else
      memory.emplace_back(std::max<size_t>(n_ele, memory.back().size() * 2L));
#endif
      it = std::prev(memory.end());
    }

    cur_head = iterator{it};
  }

public:
  simple_mem_stack() {
    clear();
  }

  /// pointers etc. are invalided on copy so there is no point in doing any work
  /// here (i.e. allocating any memory that is not used yet)
  simple_mem_stack(const simple_mem_stack&) {
    clear();
  }
  /// nothing is invalidated so this fine. Move over all the memory and keep
  /// the current head as is
  simple_mem_stack(simple_mem_stack &&o) = default;

  /// clears the object deallocating all memory
  void clear(){
    while(!marks.empty())
      marks.pop();
    memory.clear();
    memory.emplace_back(min_block_size);
    cur_head = iterator{memory.begin()};
  }

  /// returns a given number of units of memory
  T* get(size_t const n_ele) {
    if(cur_head + n_ele >= cur_head.cont->end())
      new_block(n_ele);

    T* res = &*cur_head;
    cur_head += n_ele;
    return res;
  }

  /// sets a mark in the memory that can be returned to in the future
  void set_mark(){
    marks.emplace(cur_head);
  }

  /// turn back the memory to the previous mark or the start if there is not any
  void reset_to_mark() {
    if(!marks.empty())
      cur_head = marks.top();
    else
      cur_head = iterator{memory.begin()};
  }

  /// removes the last mark
  void pop_mark() {
    if(!marks.empty())
      marks.pop();
  }

  /// turns back the memory to the start without deallocating the memory
  void reset(){
#ifdef DEBUG_SIMPLE_MEM_STACK
    clear();
#else
    while(!marks.empty())
      marks.pop();

    cur_head = iterator{memory.begin()};
#endif // #ifdef DEBUG_SIMPLE_MEM_STACK
  }

  /**
   * class used for RAII. It pops the marker and returns to the previous
   * marker when the object goes out of scope unless the current marker is not
   * the one used when the object was created
   */
  class return_memory_handler;
  friend class return_marker;
  class return_memory_handler {
    simple_mem_stack &mem_obj;
    iterator marker;

  public:
    return_memory_handler(simple_mem_stack &mem_obj, iterator marker):
      mem_obj{mem_obj}, marker{marker} { }

    ~return_memory_handler() {
      if(mem_obj.marks.empty() || mem_obj.marks.top() != marker)
        return;
      mem_obj.pop_mark();
      mem_obj.reset_to_mark();
    }
  };

  /**
   * sets a mark in the memory that can be returned to in the future. When the
   * returned object is destroyed, the marker is removed and the memory returned
   * to the previous marker.
   */
  return_memory_handler set_mark_raii(){
    marks.emplace(cur_head);
    return { *this, cur_head };
  }
};

} // namespace ghqCpp

#endif
