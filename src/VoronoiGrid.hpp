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
 * @file VoronoiGrid.hpp
 *
 * @brief Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIGRID_HPP
#define VORONOIGRID_HPP

#include "Box.hpp"

#include <ostream>
#include <vector>

class PointLocations;
class VoronoiCell;

/**
 * @brief Voronoi grid.
 */
class VoronoiGrid {
private:
  /*! @brief Bounding box containing the grid. */
  Box _box;

  /*! @brief Periodicity flags for the bounding box. */
  CoordinateVector< bool > _periodic;

  /*! @brief Cells of the grid. */
  std::vector< VoronoiCell * > _cells;

  /*! @brief Positions of the cell generators (in m). */
  std::vector< CoordinateVector<> > _generator_positions;

  /*! @brief PointLocations object used for fast neighbour searching. */
  PointLocations *_pointlocations;

public:
  VoronoiGrid(Box box, CoordinateVector< bool > periodic =
                           CoordinateVector< bool >(false),
              unsigned int numcell = 0);

  ~VoronoiGrid();

  unsigned int add_cell(CoordinateVector<> generator_position);
  void compute_grid();
  void finalize();

  double get_volume(unsigned int index) const;
  const CoordinateVector<> &get_centroid(unsigned int index) const;
  const CoordinateVector<> &get_generator(unsigned int index) const;
  CoordinateVector<> get_wall_normal(unsigned int wallindex) const;
  const std::vector< std::tuple< double, CoordinateVector<>, unsigned int > > &
  get_faces(unsigned int index) const;
  unsigned int get_index(const CoordinateVector<> &position) const;

  void print_grid(std::ostream &stream);
};

#endif // VORONOIGRID_HPP
