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

  // TODO: Make these template parameters? They could be std::array, but the
  // code using the projection currently uses std::vector and we want to
  // avoid copys.
  typedef std::vector<Coordinate> Corners;
  typedef std::vector<Coordinate> Normals;

  typedef std::array<Coordinate, dim> Images;
  typedef Images Preimages;

//private:
  const std::pair<Corners, Corners>& m_corners;
  const std::pair<Normals, Normals>& m_normals;

  /**
   * Overlap allowed for the projection to be considered valid.
   */
  const Field m_overlap = Field(0);

  /**
   * epsilon used for floating-point comparisons
   *
   * \ref epsilon(Field)
   */
  Field m_epsilon = Field(1e-4);

  /** \copydoc images() */
  std::pair<Images, Preimages> m_images;

  /** \copydoc success() */
  std::pair<std::bitset<dim>, std::bitset<dim> > m_success;

  /** \copydoc numberOfEdgeIntersections() */
  unsigned m_number_of_edge_intersections;

  /** \copydoc edgeIntersections() */
  std::array<EdgeIntersection, maxEdgeIntersections> m_edge_intersections;

  /**
   * If <code>true</code>, the forward projection was successful, that is
   * Φ(xᵢ) could be computed for all xᵢ.
   *
   * \warning Note that this only means Φ(xᵢ) lie in the plane spanned by the
   *          image simplex which is required to compute the inverse
   *          projection Φ⁻¹(yᵢ). The bitset \ref m_success should be used to
   *          check whether the projection is feasible.
   */
  bool m_projection_valid;

  /**
   * Compute forward projection Φ(xᵢ) for all xᵢ.
   */
  void doProjection();

  /**
   * Compute inverse projection Φ(yᵢ) for all yᵢ.
   *
   * \note This requires the forward projection was already computed by
   *       \ref doProjection.
   */
  void doInverseProjection();

  /**
   * Compute intersections between projected edges and edges of the image simplex.
   *
   * \note This requires the forward and inverse projections were already
   *       computed by \ref doProjection and \ref doInverseProjection.
   */
  void doEdgeIntersection();

  /**
   * Check if projection is feasible.
   *
   * Given a point <code>x</code>, its image <code>px</code> in barycentric
   * coordinates together with the signed distance along the normal at
   * <code>x</code> in the last entry of <code>px</code> and the corners and
   * normals of the image simplex given in <code>corners</code> and
   * <code>normals</code>, this method checks that the projection is feasible.
   * This means:
   *
   * <ul>
   *   <li><code>px</code> is inside the image simplex</li>
   *   <li>The signed distance given is not smaller than <code>-\ref m_overlap</code></li>
   *   <li>The signed distance along the normal at <code>px</code> is not smaller than <code>-\ref m_overlap</li>
   * </ul>
   *
   * \param x euclidean coordinate of point to project
   * \param px barycentric coordinates of projected point;
   *           last entry is distance along normal
   * \param corners corners of image simplex
   * \param normals normals of image simplex
   * \return <code>true</code> if the projection is feasible, <code>false</code> otherwise.
   */
  inline bool projectionFeasible(const Coordinate& x, const Coordinate& px, const Corners& corners, const Normals& normals) const;

public:
  // TODO: Pass corners and normals to project to avoid the reference?
  /**
   * \warning This class stores a reference to <code>corners</code> and <code>normals</code> to avoid copies.
   *          Both have to stay valid until \ref project() is called.
   *
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
   * The first d-1 values are the barycentric coordinates with respect
   * to the corners of the (pre)image, the last value is the signed
   * distance between the projected point and its (pre)image.
   *
   * \note \ref project() must be called before this method can be used.
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
   * \note \ref project() must be called before this method can be used.
   *
   * \returns pair of bitsets indicating success of (inverse) projection at
   *          corners <code>xᵢ</code> (<code>yⱼ</code>)
   */
  const std::pair<std::bitset<dim>, std::bitset<dim> >& success() const
    { return m_success; }

  /**
   * Number of edge intersections.
   *
   * \note \ref project() must be called before this method can be used.
   *
   * \ref edgeIntersections()
   */
  unsigned numberOfEdgeIntersections() const
    { return m_number_of_edge_intersections; }

  /**
   * Edge intersections.
   *
   * \warning Only the first \ref numberOfEdgeIntersections() entries are valid
   *          edge intersections.
   *
   * \seealso EdgeIntersection
   */
  const std::array<EdgeIntersection, maxEdgeIntersections>& edgeIntersections() const
    { return m_edge_intersections; }
};

} /* namespace GridGlue */
} /* namespace Dune */

#include "projectionhelper2.t.hh"

#endif
