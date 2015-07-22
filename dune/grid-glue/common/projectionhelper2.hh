#ifndef DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER2_HH
#define DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER2_HH

#include <array>
#include <bitset>
#include <tuple>
#include <utility>

namespace Dune {
namespace GridGlue {

template<typename Coordinate>
class Projection
{
public:
  constexpr static unsigned dim = Coordinate::dimension;
  constexpr static unsigned maxEdgeIntersections = dim == 3 ? 9 : 0;

  static_assert(dim == 2 || dim == 3, "Projection only implemented for dim=2 and dim=3");

  typedef typename Coordinate::field_type Field;
  typedef std::array<Coordinate, dim> Corners;
  typedef std::array<Coordinate, dim> Normals;

  const std::pair<Corners, Corners>& m_corners;
  const std::pair<Normals, Normals>& m_normals;
  const Field m_overlap = Field(0);
  const Field m_epsilon = 1e-6;
  std::pair<Corners, Corners> m_images;
  std::pair<std::bitset<dim>, std::bitset<dim> > m_success;
  unsigned m_number_of_edge_intersections;
  std::array<std::pair<Coordinate, Coordinate>, maxEdgeIntersections> m_edge_intersections;
  bool m_projection_valid;

  void doProjection();
  void doInverseProjection();
  void doEdgeIntersection();

public:
  Projection(const std::pair<Corners, Corners>& corners, const std::pair<Normals, Normals>& normals);
};

} /* namespace GridGlue */
} /* namespace Dune */

#include "projectionhelper2.t.hh"

#endif
