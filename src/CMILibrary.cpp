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
 * @file CMILibrary.cpp
 *
 * @brief CMacIonize library exposure: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CMILibrary.hpp"
#include "IonizationSimulation.hpp"
#include "SPHArrayDensityFunction.hpp"
#include "SPHArrayDensityGridWriter.hpp"

IonizationSimulation *global_ionization_simulation = nullptr;
SPHArrayDensityFunction *global_density_function = nullptr;
SPHArrayDensityGridWriter *global_density_grid_writer = nullptr;

/**
 * @brief Initialize the CMI library.
 *
 * Non-periodic version.
 *
 * @param parameter_file Parameter file to use.
 * @param num_thread Number of shared memory parallel threads to use.
 * @param unit_length_in_SI Length unit used internally (in m).
 * @param unit_mass_in_SI Mass unit used internally (in kg).
 */
void cmi_init(const char *parameter_file, const int num_thread,
              const double unit_length_in_SI, const double unit_mass_in_SI) {

  global_ionization_simulation = new IonizationSimulation(
      false, false, false, num_thread, parameter_file, nullptr, nullptr);
  global_density_function =
      new SPHArrayDensityFunction(unit_length_in_SI, unit_mass_in_SI);
  global_density_grid_writer = new SPHArrayDensityGridWriter();
}

/**
 * @brief Initialize the CMI library.
 *
 * Periodic double precision version.
 *
 * @param parameter_file Parameter file to use.
 * @param num_thread Number of shared memory parallel threads to use.
 * @param unit_length_in_SI Length unit used internally (in m).
 * @param unit_mass_in_SI Mass unit used internally (in kg).
 * @param box_anchor Coordinates of the left front lower corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 */
void cmi_init_periodic_dp(const char *parameter_file, const int num_thread,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI,
                          const double *box_anchor, const double *box_sides) {

  global_ionization_simulation = new IonizationSimulation(
      false, false, false, num_thread, parameter_file, nullptr, nullptr);
  global_density_function = new SPHArrayDensityFunction(
      unit_length_in_SI, unit_mass_in_SI, box_anchor, box_sides);
  global_density_grid_writer = new SPHArrayDensityGridWriter();
}

/**
 * @brief Initialize the CMI library.
 *
 * Periodic single precision version.
 *
 * @param parameter_file Parameter file to use.
 * @param num_thread Number of shared memory parallel threads to use.
 * @param unit_length_in_SI Length unit used internally (in m).
 * @param unit_mass_in_SI Mass unit used internally (in kg).
 * @param box_anchor Coordinates of the left front lower corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 */
void cmi_init_periodic_sp(const char *parameter_file, const int num_thread,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI, const float *box_anchor,
                          const float *box_sides) {

  global_ionization_simulation = new IonizationSimulation(
      false, false, false, num_thread, parameter_file, nullptr, nullptr);
  global_density_function = new SPHArrayDensityFunction(
      unit_length_in_SI, unit_mass_in_SI, box_anchor, box_sides);
  global_density_grid_writer = new SPHArrayDensityGridWriter();
}

/**
 * @brief Free up the memory used by the CMI library.
 */
void cmi_destroy() {
  delete global_ionization_simulation;
  delete global_density_function;
  delete global_density_grid_writer;
}

/**
 * @brief Compute the neutral fractions for the given SPH density field and
 * store them in the given array.
 *
 * Double precision version.
 *
 * @param x X coordinates (in internal units).
 * @param y Y coordinates (in internal units).
 * @param z Z coordinates (in internal units).
 * @param h Smoothing lengths (in internal units).
 * @param m Masses (in internal units).
 * @param nH Array to store the resulting neutral fractions in.
 * @param N Number of elements in each array.
 */
void cmi_compute_neutral_fraction_dp(const double *x, const double *y,
                                     const double *z, const double *h,
                                     const double *m, double *nH,
                                     const size_t N) {

  global_density_function->reset(x, y, z, h, m, N);
  global_ionization_simulation->initialize(global_density_function);
  global_density_grid_writer->reset(N, global_density_function->get_octree());
  global_ionization_simulation->run(global_density_grid_writer);
  global_density_grid_writer->fill_array(nH);
}

/**
 * @brief Compute the neutral fractions for the given SPH density field and
 * store them in the given array.
 *
 * Mixed precision version.
 *
 * @param x X coordinates (in internal units).
 * @param y Y coordinates (in internal units).
 * @param z Z coordinates (in internal units).
 * @param h Smoothing lengths (in internal units).
 * @param m Masses (in internal units).
 * @param nH Array to store the resulting neutral fractions in.
 * @param N Number of elements in each array.
 */
void cmi_compute_neutral_fraction_mp(const double *x, const double *y,
                                     const double *z, const float *h,
                                     const float *m, float *nH,
                                     const size_t N) {

  global_density_function->reset(x, y, z, h, m, N);
  global_ionization_simulation->initialize(global_density_function);
  global_density_grid_writer->reset(N, global_density_function->get_octree());
  global_ionization_simulation->run(global_density_grid_writer);
  global_density_grid_writer->fill_array(nH);
}

/**
 * @brief Compute the neutral fractions for the given SPH density field and
 * store them in the given array.
 *
 * Single precision version.
 *
 * @param x X coordinates (in internal units).
 * @param y Y coordinates (in internal units).
 * @param z Z coordinates (in internal units).
 * @param h Smoothing lengths (in internal units).
 * @param m Masses (in internal units).
 * @param nH Array to store the resulting neutral fractions in.
 * @param N Number of elements in each array.
 */
void cmi_compute_neutral_fraction_sp(const float *x, const float *y,
                                     const float *z, const float *h,
                                     const float *m, float *nH,
                                     const size_t N) {

  global_density_function->reset(x, y, z, h, m, N);
  global_ionization_simulation->initialize(global_density_function);
  global_density_grid_writer->reset(N, global_density_function->get_octree());
  global_ionization_simulation->run(global_density_grid_writer);
  global_density_grid_writer->fill_array(nH);
}
