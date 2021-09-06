#ifndef PROF_H_VAJOINT
#define PROF_H_VAJOINT
#include <string>
#include <atomic>

class profiler {
#ifdef DO_PROF_VAJOINT
  static std::atomic<bool> running_profiler;
#endif

public:
  profiler(const std::string&);

#ifdef DO_PROF_VAJOINT
  ~profiler();
#endif
};

#ifdef DO_PROF_VAJOINT
#include <gperftools/profiler.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include "psqn-Rcpp-wrapper.h"

std::atomic<bool> profiler::running_profiler(false);

profiler::~profiler(){
  ProfilerStop();
  running_profiler = false;
}
#endif

profiler::profiler(const std::string &name)
{
#ifdef DO_PROF_VAJOINT
  if(running_profiler)
    Rcpp::stop("Already running profiler...");
  running_profiler = true;

  std::stringstream ss;
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  ss << "profile-" << name << std::put_time(&tm, "-%d-%m-%Y-%H-%M-%S.log");
  Rcpp::Rcout << "Saving profile output to '" << ss.str() << "'" << std::endl;
  const std::string s = ss.str();
  ProfilerStart(s.c_str());
#endif
}

#endif
