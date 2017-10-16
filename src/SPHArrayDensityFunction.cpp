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
 * @file SPHArrayDensityFunction.cpp
 *
 * @brief SPHArrayDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "SPHArrayDensityFunction.hpp"
#include <cfloat>

/**
 * @brief Cubic spline kernel used in Gadget2.
 *
 * @param u Distance in units of the smoothing length.
 * @param h Smoothing length.
 * @return Value of the cubic spline kernel.
 */
double SPHArrayDensityFunction::cubic_spline_kernel(double u, double h) {
  const double KC1 = 2.546479089470;
  const double KC2 = 15.278874536822;
  const double KC5 = 5.092958178941;
  if (u < 1.) {
    if (u < 0.5) {
      return (KC1 + KC2 * (u - 1.) * u * u) / (h * h * h);
    } else {
      return KC5 * (1. - u) * (1. - u) * (1. - u) / (h * h * h);
    }
  } else {
    // the cubic spline kernel has compact support
    return 0.;
  }
}

/**
 * @brief Constructor.
 *
 * Non-periodic version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 */
SPHArrayDensityFunction::SPHArrayDensityFunction(const double unit_length_in_SI,
                                                 const double unit_mass_in_SI)
    : _unit_length_in_SI(unit_length_in_SI), _unit_mass_in_SI(unit_mass_in_SI),
      _is_periodic(false), _octree(nullptr) {}

/**
 * @brief Constructor.
 *
 * Periodic double precision version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param box_anchor Coordinates of the left front bottom corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 */
SPHArrayDensityFunction::SPHArrayDensityFunction(const double unit_length_in_SI,
                                                 const double unit_mass_in_SI,
                                                 const double *box_anchor,
                                                 const double *box_sides)
    : _unit_length_in_SI(unit_length_in_SI), _unit_mass_in_SI(unit_mass_in_SI),
      _is_periodic(true), _octree(nullptr) {

  _box.get_anchor()[0] = box_anchor[0] * _unit_length_in_SI;
  _box.get_anchor()[1] = box_anchor[1] * _unit_length_in_SI;
  _box.get_anchor()[2] = box_anchor[2] * _unit_length_in_SI;
  _box.get_sides()[0] = box_sides[0] * _unit_length_in_SI;
  _box.get_sides()[1] = box_sides[1] * _unit_length_in_SI;
  _box.get_sides()[2] = box_sides[2] * _unit_length_in_SI;
}

/**
 * @brief Constructor.
 *
 * Periodic single precision version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param box_anchor Coordinates of the left front bottom corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 */
SPHArrayDensityFunction::SPHArrayDensityFunction(const double unit_length_in_SI,
                                                 const double unit_mass_in_SI,
                                                 const float *box_anchor,
                                                 const float *box_sides)
    : _unit_length_in_SI(unit_length_in_SI), _unit_mass_in_SI(unit_mass_in_SI),
      _is_periodic(true), _octree(nullptr) {

  _box.get_anchor()[0] = box_anchor[0] * _unit_length_in_SI;
  _box.get_anchor()[1] = box_anchor[1] * _unit_length_in_SI;
  _box.get_anchor()[2] = box_anchor[2] * _unit_length_in_SI;
  _box.get_sides()[0] = box_sides[0] * _unit_length_in_SI;
  _box.get_sides()[1] = box_sides[1] * _unit_length_in_SI;
  _box.get_sides()[2] = box_sides[2] * _unit_length_in_SI;
}

/**
 * @brief Destructor.
 *
 * Frees up memory used by the internal Octree.
 */
SPHArrayDensityFunction::~SPHArrayDensityFunction() { delete _octree; }

/**
 * @brief Reset the internal data values.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayDensityFunction::reset(const double *x, const double *y,
                                    const double *z, const double *h,
                                    const double *m, const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Reset the internal data values.
 *
 * This version uses single precision smoothing length and masses.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayDensityFunction::reset(const double *x, const double *y,
                                    const double *z, const float *h,
                                    const float *m, const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Reset the internal data values.
 *
 * This version uses single precision coordinates, smoothing lengths and masses.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayDensityFunction::reset(const float *x, const float *y,
                                    const float *z, const float *h,
                                    const float *m, const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Get a pointer to the internal Octree.
 *
 * @return Pointer to the internal Octree.
 */
Octree *SPHArrayDensityFunction::get_octree() { return _octree; }

/**
 * @brief Initialize the internal Octree.
 */
void SPHArrayDensityFunction::initialize() {
  _octree = new Octree(_positions, _box, _is_periodic);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);
}


/**
 * @brief Function that gives the 3D integral of the kernel of
 * a particle for a given vertex of a cell face.
 *
 * @param phi Azimuthal angle of the vertex.
 * @param r0 Distance from the particle to the face of the cell.
 * @param R_0 Distance from the orthogonal projection of the particle
 * onto the face of the cell to a side of the face (containing the vertex).
 * @param h The kernel smoothing length of the particle.
 * @return The integral of the kernel for the given vertex.
 */

double SPHArrayDensityFunction::full_integral(double phi, double r0,
                                                   double R_0, double h) {

  double B1, B2, B3, mu, a, logs, u;
  double full_int;
  double a2, cosp, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3, tanp;
  double r2, R, linedist2, phi1, phi2, h2;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double D2, D3;

  B1 = 0.0;
  B2 = 0.0;
  B3 = 0.0;

  if (r0 == 0.0)
    return 0.0;
  if (R_0 == 0.0)
    return 0.0;
  if (phi == 0.0)
    return 0.0;

  h2 = h * h;
  r03 = r0 * r0 * r0;
  r0h2 = r0 / h * r0 / h;
  r0h3 = r0h2 * r0 / h;
  r0h_2 = h / r0 * h / r0;
  r0h_3 = r0h_2 * h / r0;

  // Setting up the B1, B2, B3 constants of integration.

  if (r0 >= 2.0 * h) {
    B3 = h2 * h / 4.;
  } else if (r0 > h) {
    B3 = r03 / 4. * (-4. / 3. + (r0 / h) - 0.3 * r0h2 + 1. / 30. * r0h3 -
                     1. / 15. * r0h_3 + 8. / 5. * r0h_2);
    B2 = r03 / 4. * (-4. / 3. + (r0 / h) - 0.3 * r0h2 + 1. / 30. * r0h3 -
                     1. / 15. * r0h_3);
  } else {
    B3 = r03 / 4. * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3 + 7. / 5. * r0h_2);
    B2 = r03 / 4. * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3 - 1. / 5. * r0h_2);
    B1 = r03 / 4. * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3);
  }

  a = R_0 / r0;
  a2 = a * a;

  linedist2 = r0 * r0 + R_0 * R_0;
  R = R_0 / cos(phi);
  r2 = r0 * r0 + R * R;

  full_int = 0.0;
  D2 = 0.0;
  D3 = 0.0;

  if (linedist2 <= h2) {
    ////// phi1 business /////
    phi1 = acos(R_0 / sqrt(h * h - r0 * r0));

    cosp = cos(phi1);
    cosp2 = cosp * cosp;
    mu = cosp / a / sqrt(1. + cosp2 / a2);

    tanp = tan(phi1);

    I0 = phi1;
    I_2 = phi1 + a2 * tanp;
    I_4 = phi1 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sin(phi1) * sqrt(1 - mu * mu);
    logs = log(1 + u) - log(1 - u);
    I1 = atan(u / a);

    I_1 = a / 2. * logs + I1;
    I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    I_5 =
        I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

    D2 = -1. / 6. * I_2 + 0.25 * (r0 / h) * I_3 - 0.15 * r0h2 * I_4 +
         1. / 30. * r0h3 * I_5 - 1. / 60. * r0h_3 * I1 + (B1 - B2) / r03 * I0;

    ////// phi2 business /////
    phi2 = acos(R_0 / sqrt(4.0 * h * h - r0 * r0));

    cosp = cos(phi2);
    cosp2 = cosp * cosp;
    mu = cosp / a / sqrt(1. + cosp2 / a2);

    tanp = tan(phi2);

    I0 = phi2;
    I_2 = phi2 + a2 * tanp;
    I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sin(phi2) * sqrt(1 - mu * mu);
    logs = log(1 + u) - log(1 - u);
    I1 = atan(u / a);

    I_1 = a / 2. * logs + I1;
    I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    I_5 =
        I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

    D3 = 1. / 3. * I_2 - 0.25 * (r0 / h) * I_3 + 3. / 40. * r0h2 * I_4 -
         1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 +
         D2;
  } else if (linedist2 <= 4.0 * h2) {
    ////// phi2 business /////
    phi2 = acos(R_0 / sqrt(4.0 * h * h - r0 * r0));

    cosp = cos(phi2);
    cosp2 = cosp * cosp;
    mu = cosp / a / sqrt(1. + cosp2 / a2);

    tanp = tan(phi2);

    I0 = phi2;
    I_2 = phi2 + a2 * tanp;
    I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sin(phi2) * sqrt(1 - mu * mu);
    logs = log(1 + u) - log(1 - u);
    I1 = atan(u / a);

    I_1 = a / 2. * logs + I1;
    I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    I_5 =
        I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

    D3 = 1. / 3. * I_2 - 0.25 * (r0 / h) * I_3 + 3. / 40. * r0h2 * I_4 -
         1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 +
         D2;
  }

  //////////////////////////////////
  // Calculating I_n expressions. //
  //////////////////////////////////

  cosp = cos(phi);
  cosp2 = cosp * cosp;
  mu = cosp / a / sqrt(1. + cosp2 / a2);

  tanp = tan(phi);

  I0 = phi;
  I_2 = phi + a2 * tanp;
  I_4 = phi + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

  u = sin(phi) * sqrt(1 - mu * mu);
  logs = log(1 + u) - log(1 - u);
  I1 = atan(u / a);

  I_1 = a / 2. * logs + I1;
  I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
  I_5 = I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

  // Calculating the integral expression.

  if (r2 < h2) {
    full_int = r0h3 / M_PI * (1. / 6. * I_2 - 3. / 40. * r0h2 * I_4 +
                              1. / 40. * r0h3 * I_5 + B1 / r03 * I0);
  } else if (r2 < 4.0 * h2) {
    full_int = r0h3 / M_PI *
               (0.25 * (4. / 3. * I_2 - (r0 / h) * I_3 + 0.3 * r0h2 * I_4 -
                        1. / 30. * r0h3 * I_5 + 1. / 15. * r0h_3 * I1) +
                B2 / r03 * I0 + D2);
  } else {
    full_int = r0h3 / M_PI * (-0.25 * r0h_3 * I1 + B3 / r03 * I0 + D3);
  }

  return full_int;
}

/**
 * @brief Function that calculates the mass contribution of a particle
 * towards the total mass of a cell.
 *
 * @param cell Geometrical information about the cell.
 * @param particle The particle position.
 * @param h The kernel smoothing length of the particle.
 * @return The mass contribution of the particle to the cell.
 */

double SPHArrayDensityFunction::mass_contribution(
    const Cell &cell, const CoordinateVector<> particle, const double h) {
  double M, Msum;

  Msum = 0.;
  M = 0.;

  std::vector< Face > face_vector = cell.get_faces();

  // Loop over each face of a cell.
  for (unsigned int i = 0; i < face_vector.size(); i++) {

    CoordinateVector<> vert_position1;
    CoordinateVector<> projected_particle;
    double r0 = 0.;
    double ar0 = 0.;
    double s2 = 0.;

    // Loop over the vertices of each face.
    for (Face::Vertices j = face_vector[i].first_vertex();
         j != face_vector[i].last_vertex(); ++j) {
      if (j == face_vector[i].first_vertex()) { // Calculating the distance from
                                                // particle to each face of a
                                                // cell
        // http://mathinsight.org/distance_point_plane
        // http://mathinsight.org/forming_planes
        Face::Vertices j_twin = j;
        vert_position1 = j_twin.get_position();
        CoordinateVector<> vert_position2 = (++j_twin).get_position();
        CoordinateVector<> vert_position3 = (++j_twin).get_position();

        const double A = (vert_position2[1] - vert_position1[1]) *
                             (vert_position3[2] - vert_position1[2]) -
                         (vert_position3[1] - vert_position1[1]) *
                             (vert_position2[2] - vert_position1[2]);
        const double B = (vert_position2[2] - vert_position1[2]) *
                             (vert_position3[0] - vert_position1[0]) -
                         (vert_position3[2] - vert_position1[2]) *
                             (vert_position2[0] - vert_position1[0]);
        const double C = (vert_position2[0] - vert_position1[0]) *
                             (vert_position3[1] - vert_position1[1]) -
                         (vert_position3[0] - vert_position1[0]) *
                             (vert_position2[1] - vert_position1[1]);
        const double D = -A * vert_position1[0] - B * vert_position1[1] -
                         C * vert_position1[2];

        const double norm = sqrt(A * A + B * B + C * C);
        r0 = (A * particle[0] + B * particle[1] + C * particle[2] + D) / norm;
        ar0 = fabs(r0);

        // Calculate of the orthogonal projection of the particle position onto
        // the face.
        projected_particle[0] = particle[0] - r0 * A / norm;
        projected_particle[1] = particle[1] - r0 * B / norm;
        projected_particle[2] = particle[2] - r0 * C / norm;

        // s2 contains information about the orientation of the face vertices.
        s2 = vert_position1[0] * (vert_position2[1] * vert_position3[2] -
                                  vert_position2[2] * vert_position3[1]) +
             vert_position1[1] * (vert_position2[2] * vert_position3[0] -
                                  vert_position2[0] * vert_position3[2]) +
             vert_position1[2] * (vert_position2[0] * vert_position3[1] -
                                  vert_position2[1] * vert_position3[0]);
      }

      Face::Vertices j_twin = j;
      const CoordinateVector<> vert_position2 = j_twin.get_position();
      CoordinateVector<> vert_position3;

      if (++j_twin == face_vector[i].last_vertex()) {
        vert_position3 = vert_position1;
      } else {
        vert_position3 = j_twin.get_position();
      }

      const double r23 = (vert_position2 - vert_position3).norm();
      const double r12 = (projected_particle - vert_position2).norm();
      const double r13 = (projected_particle - vert_position3).norm();
      const double cosa = ((vert_position3[0] - vert_position2[0]) *
                               (projected_particle[0] - vert_position2[0]) +
                           (vert_position3[1] - vert_position2[1]) *
                               (projected_particle[1] - vert_position2[1]) +
                           (vert_position3[2] - vert_position2[2]) *
                               (projected_particle[2] - vert_position2[2])) /
                          r12 / r23;

      double R_0 = 0.;
      double phi1 = 0.;
      double phi2 = 0.;

      if (fabs(cosa) < 1.0) {
        R_0 = r12 * sqrt(1 - cosa * cosa);
      } else {
        if (fabs(cosa) - 1.0 < 0.00001) {
          R_0 = 0.0;
        } else {
          printf("Error: cosa > 1: %g\n", cosa);
        }
      }

      const double s1 =
          projected_particle[0] * (vert_position2[1] * vert_position3[2] -
                                   vert_position2[2] * vert_position3[1]) +
          projected_particle[1] * (vert_position2[2] * vert_position3[0] -
                                   vert_position2[0] * vert_position3[2]) +
          projected_particle[2] * (vert_position2[0] * vert_position3[1] -
                                   vert_position2[1] * vert_position3[0]);

      if (R_0 < r12) {
        phi1 = acos(R_0 / r12);
      } else {
        if ((R_0 - r12) / h < 0.00001) {
          phi1 = 0.0;
        } else {
          printf("Error: R0 > r12: %g\n", R_0 - r12);
        }
      }
      if (R_0 < r13) {
        phi2 = acos(R_0 / r13);
      } else {
        if ((R_0 - r13) / h < 0.00001) {
          phi2 = 0.0;
        } else {
          printf("Error: R0 > r13: %g\n", R_0 - r13);
        }
      }

      // Find out if the vertex integral will contribute positively or
      // negatively to the cell mass.
      if (s1 * s2 * r0 <= 0) {
        M = -1.;
      } else {
        M = 1.;
      }

      // Calculate the vertex integral.
      if ((r12 * sin(phi1) >= r23) || (r13 * sin(phi2) >= r23)) {
        if (phi1 >= phi2) {
          M = M * (full_integral(phi1, ar0, R_0, h) -
                   full_integral(phi2, ar0, R_0, h));
        } else {
          M = M * (full_integral(phi2, ar0, R_0, h) -
                   full_integral(phi1, ar0, R_0, h));
        }
      } else {
        M = M * (full_integral(phi1, ar0, R_0, h) +
                 full_integral(phi2, ar0, R_0, h));
      }
      Msum = Msum + M;
    }
  }

  return Msum;
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues SPHArrayDensityFunction::operator()(const Cell &cell) const {

  DensityValues values;
  bool new_density_algorithm = true;

  if(new_density_algorithm) {
    CoordinateVector<> position = cell.get_cell_midpoint();

    // Find the vertex that is furthest away from the cell midpoint.
    std::vector< Face > face_vector = cell.get_faces();
    double radius = 0.0;
    for (unsigned int i = 0; i < face_vector.size(); i++) {
      for (Face::Vertices j = face_vector[i].first_vertex();
           j != face_vector[i].last_vertex(); ++j) {
        double distance = j.get_position().norm();
        if (distance > radius)
          radius = distance;
      }
    }

    // Find the neighbours that are contained inside of a sphere of centre the
    // cell midpoint
    // and radius given by the distance to the furthest vertex.
    std::vector< unsigned int > ngbs =
        _octree->get_ngbs_sphere(position, radius);
    const unsigned int numngbs = ngbs.size();

    double density = 0.;

    // Loop over all the neighbouring particles and calculate their mass
    // contributions.
    for (unsigned int i = 0; i < numngbs; i++) {
      const unsigned int index = ngbs[i];
      const double h = _smoothing_lengths[index] / 2.0;
      const CoordinateVector<> particle = _positions[index];
      density += mass_contribution(cell, particle, h) * _masses[index];
    }

    // Divide the cell mass by the cell volume to get density.
    density = density / cell.get_volume();

    // convert density to particle density (assuming hydrogen only)
    values.set_number_density(density / 1.6737236e-27);
    // TODO: other quantities
    // temporary values
    values.set_temperature(8000.);
    values.set_ionic_fraction(ION_H_n, 1.e-6);
    values.set_ionic_fraction(ION_He_n, 1.e-6);

  } else {

    const CoordinateVector<> position = cell.get_cell_midpoint();

    double density = 0.;
    const std::vector< unsigned int > ngbs = _octree->get_ngbs(position);
    const unsigned int numngbs = ngbs.size();
    for (unsigned int i = 0; i < numngbs; ++i) {
      const unsigned int index = ngbs[i];
      double r;
      if (!_box.get_sides().x()) {
        r = (position - _positions[index]).norm();
      } else {
        r = _box.periodic_distance(position, _positions[index]).norm();
      }
      const double h = _smoothing_lengths[index];
      const double u = r / h;
      const double m = _masses[index];
      const double splineval = m * cubic_spline_kernel(u, h);
      density += splineval;
    }

    values.set_number_density(density / 1.6737236e-27);
    values.set_temperature(8000.);
    values.set_ionic_fraction(ION_H_n, 1.e-6);
    values.set_ionic_fraction(ION_He_n, 1.e-6);
  }
  return values;
}
