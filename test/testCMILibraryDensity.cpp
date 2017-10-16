/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file testCMILibrary.cpp
 *
 * @brief Unit test for the CMILibrary density mapping.
 *
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#include "CMILibrary.hpp"
#include <fstream>
#include <vector>

/**
 * @brief Unit test for the CMILibrary density mapping.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // set up a test
  const double pc = 3.086e16;
  std::vector< double > x(2, 0.);
  std::vector< double > y(2, 0.);
  std::vector< double > z(2, 0.);
  std::vector< double > h(2, 0.);
  std::vector< double > m(2, 0.);
  const double box_anchor[3] = {-5. * pc, -5. * pc, -5. * pc};
  const double box_sides[3] = {10. * pc, 10. * pc, 10. * pc};

  x[0] = 1 * pc;
  y[0] = 1 * pc;
  z[0] = -1 * pc;
  h[0] = 5 * pc;
  m[0] = 1.e33;

  x[1] = -1 * pc;
  y[1] = -1 * pc;
  z[1] = 0 * pc;
  h[1] = 4 * pc;
  m[1] = 1.e33;

  // initialize the library
  cmi_init_periodic_dp("test_CMI_library_density.param", 1, 1., 1., box_anchor,
                       box_sides);

  // run the simulation
  std::vector< double > nH(2, 0.);
  cmi_compute_neutral_fraction_dp(x.data(), y.data(), z.data(), h.data(),
                                  m.data(), nH.data(), 2);

  // write an output file for visual checking
  std::ofstream ofile("test_CMI_library.txt");
  for (unsigned int i = 0; i < 2; ++i) {
    ofile << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << nH[i] << "\n";
  }
  ofile.close();

  // clean up the library
  cmi_destroy();

  return 0;
}
