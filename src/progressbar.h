
#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <cstdio>
#include <iostream>

namespace {
  std::string makeBar(std::string b) {
    std::string bar;
    for (int i = 0; i < 100; ++i) {
      bar += b;
    }
    return bar;
  }
  std::string null_bar;
  std::string fill_bar;
}  // namespace

//Header only, it is implemented in the progressbar.cpp file
bool progressBar(int i, int n, std::string fill, std::string null);
#endif
