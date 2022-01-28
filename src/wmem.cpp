#include "wmem.h"
#include <vector>

namespace wmem {
std::vector<ghqCpp::simple_mem_stack<double> > mem_stacks =
  std::vector<ghqCpp::simple_mem_stack<double> >(1);
std::vector<ghqCpp::simple_mem_stack<cfaad::Number> > mem_stacks_Num =
  std::vector<ghqCpp::simple_mem_stack<cfaad::Number> >(1);

void setup_working_memory(const size_t n_threads){
  mem_stacks_Num.resize(n_threads);
  mem_stacks.resize(n_threads);
}

void rewind(const size_t idx){
  mem_stacks_Num[idx].reset();
  mem_stacks[idx].reset();
}

void rewind_to_mark(const size_t idx){
  mem_stacks_Num[idx].reset_to_mark();
  mem_stacks[idx].reset_to_mark();
}

void set_mark(const size_t idx){
  mem_stacks_Num[idx].set_mark();
  mem_stacks[idx].set_mark();
}

void rewind_all(){
  for(auto &x : mem_stacks_Num) x.reset();
  for(auto &x : mem_stacks) x.reset();
}

void clear_all(){
  for(auto &x : mem_stacks_Num) x.clear();
  for(auto &x : mem_stacks) x.clear();
}

double * get_double_mem(const size_t n){
  return mem_stacks[get_thread_num()].get(n);
}

cfaad::Number * get_Number_mem(const size_t n){
  return mem_stacks_Num[get_thread_num()].get(n);
}

ghqCpp::simple_mem_stack<double> &mem_stack(const size_t idx){
  return mem_stacks[idx];
}

} // namespace wmem
