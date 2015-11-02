#ifndef DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER2_HH
#define DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER2_HH

#include <array>
#include <bitset>
#include <tuple>
#include <utility>
#include <vector>

namespace Dune {
namespace GridGlue {

/**
 * Projection of a line (triangle) on another line (triangle).
 */
template<typename Coordinate>
class Projection
{
public:
  /**
   * Intersection between two edges of a triangle.
   */
  struct EdgeIntersection
  {
    /**
     * Edge numbers in image and preimage triangle
     */
    std::array<unsigned, 2> edge;

    /**
     * Local coordinate of intersection point in barycentric coordinates with
     * respect to image and preimage triangle.
     */
    std::array<Coordinate, 2> local;
  };

  constexpr static unsigned dim = Coordinate::dimension;
  constexpr static unsigned maxEdgeIntersections = dim == 3 ? 9 : 0;

  static_assert(dim == 2 || dim == 3, "Projection only implemented for dim=2 or dim=3");

  typedef typename Coordinate::field_type Field;
  typedef std::vector<Coordinate> Corners;
  typedef std::vector<Coordinate> Normals;
  typedef std::array<Coordinate, dim> Images;
  typedef Images Preimages;

//private:
  const std::pair<Corners, Corners>& m_corners;
  const std::pair<Normals, Normals>& m_normals;
  const Field m_overlap = Field(0);

  /**
   * epsilon used for floating-point comparisons
   *
   * \ref epsilon(Field)
   */
  Field m_epsilon = 1e-4;

  /** \copydoc images() */
  std::pair<Images, Preimages> m_images;

  /** \copydoc success() */
  std::pair<std::bitset<dim>, std::bitset<dim> > m_success;

  /** \copydoc numberOfEdgeIntersections() */
  unsigned m_number_of_edge_intersections;

  /** \copydoc edgeIntersections() */
  std::array<EdgeIntersection, maxEdgeIntersections> m_edge_intersections;

  bool m_projection_valid;

  void doProjection();
  void doInverseProjection();
  void doEdgeIntersection();

  bool projectionFeasible(const Coordinate& x, const Coordinate& nx, const Coordinate& y, const Coordinate& ny) const;

public:
  /**
   * \param corners euclidean coordinates of corners of preimage and image
   * \param normals normals at corners of preimage and image
   */
  Projection(const std::pair<Corners, Corners>& corners, const std::pair<Normals, Normals>& normals, const Field overlap = Field(0));

  /**
   * Set epsilon used for floating-point comparisons.
   *
   * \param epsilon new epsilon used for floating-point comaprisons
   */
  void epsilon(const Field epsilon);

  /**
   * Do the actual projection.
   */
  void project();

  /**
   * Images and preimages of corners.
   *
   * Returns a pair of arrays. The first array contains the images
   * <code>Φ(xᵢ)</code> of the corners <code>xᵢ</code>. The second
   * array contains the preimages <code>Φ⁻¹(yⱼ)</code> of the
   * corners <code>yⱼ</code>.
   *
   * All values are barycentric coordinates with respect to the corners of the (pre)image.
   *
   * \returns pair of arrays giving <code>((Φ(xᵢ))ᵢ, (Φ⁻¹(yⱼ))ⱼ)</code> in barycentric coordinates
   *
   * \ref success()
   */
  const std::pair<Images, Preimages>& images() const
    { return m_images; }

  /**
   * Indicate whether projection (inverse projection) is valid for each corner or not.
   *
   * Returns a pair of bitsets. The first bitset indicates if the projection
   * <code>Φ(xᵢ)</code> is valid for each corner <code>xᵢ</code>, that is
   * that <code>Φ(xᵢ)</code> could be computed and lies in the image simplex.
   * The second bitset indicates the same for the inverse projection
   * <code>Φ⁻¹(yⱼ)</code> for the corners <code>yⱼ</code>.
   *
   * \returns pair of bitsets indicating success of (inverse) projection at
   *          corners <code>xᵢ</code> (<code>yⱼ</code>)
   */
  const std::pair<std::bitset<dim>, std::bitset<dim> >& success() const
    { return m_success; }

  /**
   * Number of edge intersections.
   *
   * \ref edgeIntersections()
   */
  unsigned numberOfEdgeIntersections() const
    { return m_number_of_edge_intersections; }

  /**
   * Edge intersections.
   *
   * \seealso numberOfEdgeIntersections()
   * \seealso EdgeIntersection
   */
  const std::array<EdgeIntersection, maxEdgeIntersections>& edgeIntersections() const
    { return m_edge_intersections; }
};

} /* namespace GridGlue */
} /* namespace Dune */

#include "projectionhelper2.t.hh"

#endif
