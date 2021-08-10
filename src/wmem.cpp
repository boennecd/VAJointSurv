#include "wmem.h"
#include <vector>

namespace wmem {
using double_block = cfaad::blocklist<double       , max_ele>;
using Number_block = cfaad::blocklist<cfaad::Number, max_ele>;
std::vector<double_block> double_mem = std::vector<double_block>(1);
std::vector<Number_block> Number_mem = std::vector<Number_block>(1);

void setup_working_memory(const size_t n_threads){
  double_mem.resize(n_threads);
  Number_mem.resize(n_threads);
}

void rewind(const size_t idx){
  double_mem[idx].rewind();
  Number_mem[idx].rewind();
}

void clear_all(){
  for(auto &x : double_mem) x.clear();
  for(auto &x : Number_mem) x.clear();
}

double * get_double_mem(const size_t n){
  return double_mem[get_thread_num()].emplace_back_multi(n);
}

cfaad::Number * get_Number_mem(const size_t n){
  return Number_mem[get_thread_num()].emplace_back_multi(n);
}

} // namespace wmem
