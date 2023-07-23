//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/20.
// modified by Davide Cenzato on 7/23/23
//

#include <iostream>
#include <chrono> 
#include "load-index.h"

std::size_t locateQuery(const auto &t_idx, const auto &pattern)
{
  //std::size_t total_occs = 0;
  //auto occs = t_idx.first->Locate(pattern);

  return t_idx.first->Locate(pattern).size();
}

std::size_t countQuery(const auto &t_idx, const auto &pattern)
{
  //std::size_t total_occs = 0;
  //auto occs = t_idx.first->Count(pattern);

  return t_idx.first->Count(pattern);
}

int main(int argc, char *argv[]) {

  std::string data_dir, patt_file, data_name;

  if( argc < 4 )
  {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }
  else
  {
    data_dir  = std::string(argv[1]);
    data_name = std::string(argv[2]);
    patt_file = std::string(argv[3]);
  }

  // Index Loading
  std::cout << "Loading sr-index\n";
  sdsl::cache_config config(true, data_dir, data_name);

  Factory factory(config, 64, true);

  //auto index = factory.make(Factory::Config{Factory::IndexEnum::RIndexSampled, 4});
  auto index = factory.make(Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedAreas, 64});
  //auto index = factory.make(Factory::Config{Factory::IndexEnum::RIndex});

  // Query patterns
  {
    std::ifstream pattern_file(patt_file, std::ios_base::binary);
    if (!pattern_file) {
      std::cerr << "ERROR: Failed to open patterns file! (" << patt_file << ")" << std::endl;
      return 2;
    }

    std::string buf;
    std::size_t occ_tot = 0, no_patt = 0;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(pattern_file, buf))
    {
      no_patt++;
      //if (buf.empty()) 
      if( no_patt%2 != 0 )
        continue;

      //std::cout << "pat: " << buf << "\n";
      //std::cout << "Locate nocc: " << locateQuery(index,buf) << "\n";
      //std::cout << "Count nocc: "  << countQuery(index,buf) << "\n";
      //patterns.emplace_back(buf);
      occ_tot += locateQuery(index,buf);
    }
    auto end = std::chrono::high_resolution_clock::now();
    pattern_file.close();

    uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    no_patt /= 2;
    std::cout << "number of patterns n = " << no_patt << std::endl;
    //cout << "pattern length m = " << m << endl;
    std::cout << "total number of occurrences  occ_t = " << occ_tot << std::endl;

    std::cout << "Total time : " << search << " milliseconds" << std::endl;
    std::cout << "Search time : " << (double)search/no_patt << " milliseconds/pattern (total: " << no_patt << " patterns)" << std::endl;
    std::cout << "Search time : " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << std::endl;
  }

  return 0;
}
