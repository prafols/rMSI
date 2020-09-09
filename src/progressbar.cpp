#include "progressbar.h"


bool progressBar(int i, int n, std::string fill, std::string null)
{
  if (null_bar.size() == 0)
    null_bar = makeBar(null);
  if (fill_bar.size() == 0)
    fill_bar = makeBar(fill);
  if (i == n) { i--; } //Ensure i is below n at 100%
  if (i < n) {
    int progress = ((i+1)*100)/n;
    printf("\r%d%%|", progress);
    std::cout << fill_bar.substr(0, fill.size()*progress)
              << null_bar.substr(0, null.size()*(100 - progress));
    printf("| %d/%d", i+1, n);
    progress == 100 ? std::cout << std::endl : std::cout << std::flush;
    return true;
  }
  null_bar = "";
  fill_bar = "";
  printf("\n"); //End with an empty line
  return false;
}