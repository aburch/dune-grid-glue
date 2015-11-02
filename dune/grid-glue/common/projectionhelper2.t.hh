#include <dune/common/fmatrix.hh>
#include <limits>

namespace Dune {
namespace GridGlue {

namespace ProjectionImplementation {

/**
 * Return corner coordinates of a simplex.
 *
 * Given the number <code>c</code> of a corner, this function returns
 * the coordinate of the <code>c</code>th corner of the standard
 * simplex with the same dimension as <code>Coordinate</code>.
 *
 * \param c corner number
 * \returns coordinate of <code>c</code>th corner of the standard simplex
 */
template<typename Coordinate, typename Field>
inline Coordinate
corner(unsigned c)
{
  Coordinate x(Field(0));
  if (c == 0)
    return x;
  x[c-1] = Field(1);
  return x;
}

/**
 * Translate edge to corner numbers.
 *
 * Given the number <code>edge</code> of an edge of a triangle, this
 * function returns the number of the corners belonging to it.
 *
 * \param edge edge number of a triangle
 * \return numbers of corners of the edge
 */
inline std::pair<unsigned, unsigned>
edgeToCorners(unsigned edge)
{
  switch(edge) {
    case 0: return {0, 1};
    case 1: return {0, 2};
    case 2: return {1, 2};
  }
  DUNE_THROW(Dune::Exception, "Unexpected edge number.");
}

/**
 * Convert barycentric coordinates to euclidian coordinates.
 *
 * This function converts barycentric coordinates <code>x</code> with respect
 * to the triangle with corners <code>corners</code> to euclidian coordinates.
 * For the result <code>y</code> the following equation holds:
 * <code>yᵢ = (cornersᵢ₊₁ - corners₀) xᵢ</code>
 *
 * Note that this can also be for linear interpolation of normals given on the corners.
 *
 * \param x barycentric coordinates
 * \param corners coordinates or normals at the corners
 * \return euclidian coordinates or interpolated normal
 */
template<typename Coordinate, typename Corners>
inline typename Corners::value_type
interpolate(const Coordinate& x, const Corners& corners)
{
  auto y = corners[0];
  for (unsigned i = 0; i < corners.size() - 1; ++i)
    y.axpy(x[i], corners[i+1] - corners[0]);
  return y;
}

/**
 * Check if the point <code>x</code> is inside the standard simplex.
 *
 * This functions checks if the point <code>x</code> is in the inside
 * (or on the boundary) of the standard simplex, that is
 * <code>xᵢ ≥ 0</code> and <code>∑ xᵢ ≤ 1</code>.
 *
 * \param x coordinates of point to check
 * \param epsilon tolerance used for floating-point comparisions
 * \return <code>true</code> if <code>x</code> is inside, <code>false</code> otherwise
 */
template<typename Coordinate, typename Field>
inline bool
inside(const Coordinate& x, const Field& epsilon)
{
  const unsigned dim = Coordinate::dimension;
  Field sum(0);
  for (unsigned i = 0; i < dim-1; ++i) {
    if (x[i] < -epsilon || x[i] > Field(1) + epsilon)
      return false;
    sum += x[i];
  }
  if (sum > Field(1) + epsilon)
    return false;
  return true;
}

#if 0
template<typename Local, typename Coordinate>
inline Coordinate::value_type
distanceAlongNormal(const Local& xlocal, const std::array<Coordinate, Coordinate::dimension>& xnormals, const Coordinate& x, const Coordinate& y)
{
  const auto& n = interpolate(xlocal, xnormals);
  return (y-x)*n;
}
#endif

} /* namespace ProjectionImplementation */

template<typename Coordinate>
Projection<Coordinate>
::Projection(const std::pair<Corners, Corners>& corners, const std::pair<Normals, Normals>& normals, const Field overlap)
  : m_corners(corners)
  , m_normals(normals)
  , m_overlap(overlap)
{
  /* Nothing. */
}

template<typename Coordinate>
void
Projection<Coordinate>
::epsilon(const Field epsilon)
{
  m_epsilon = epsilon;
}

template<typename Coordinate>
void
Projection<Coordinate>
::doProjection()
{
  using namespace ProjectionImplementation;
  typedef Dune::FieldMatrix<Field, dim, dim> Matrix;
  Matrix m;

  const auto& origin = m_corners.first;
  const auto& normals = m_normals.first;
  const auto& target = m_corners.second;
  auto& images = m_images.first;
  auto& success = m_success.first;

  std::array<Coordinate, dim-1> directions;
  std::array<Field, dim-1> scales;
  for (unsigned i = 0; i < dim-1; ++i) {
    directions[i] = target[i+1] - target[0];
    scales[i] = directions[i].infinity_norm();
    directions[i] /= scales[i];
  }

  for (unsigned i = 0; i < dim-1; ++i) {
    for (unsigned j = 0; j < dim; ++j) {
      m[j][i] = directions[i][j];
    }
  }

  m_projection_valid = true;
  success.reset();

  /* Project xᵢ */
  for (unsigned i = 0; i < origin.size(); ++i) {
    //std::cout << "do projection of " << origin[i] << std::endl;
    for (unsigned j = 0; j < dim; ++j)
      m[j][dim-1] = normals[i][j];

    const Coordinate rhs = origin[i] - target[0];
    try {
      //std::cout << "m = " << m << std::endl;
      auto& y = images[i];
      m.solve(y, rhs);
      for (unsigned j = 0; j < dim-1; ++j)
        y[j] /= scales[j];
      y[dim-1] *= Field(-1);

      bool success_i = inside(y, m_epsilon) && y[dim-1] > -m_overlap-m_epsilon;
      success.set(i, success_i);
    }
    catch (const Dune::FMatrixError&) {
      success.set(i, false);
      m_projection_valid = false;
    }
  }

#if 0
  std::cout << "-------------------------------------------" << std::endl
            << "projection: success = " << success << std::endl
            << "  origin: " << std::endl;
  for (const auto& x : origin)
    std::cout << "    " << x << std::endl;
  std::cout << "  target:" << std::endl;
  for (const auto& x : target)
    std::cout << "    " << x << std::endl;
  std::cout << std::endl;
#endif
}

template<typename Coordinate>
void
Projection<Coordinate>
::doInverseProjection()
{
  using namespace ProjectionImplementation;
  using std::get;
  typedef Dune::FieldMatrix<Field, dim-1, dim-1> Matrix;
  typedef Dune::FieldVector<Field, dim-1> Vector;

  if (!m_projection_valid) {
    m_success.second.reset();
    return;
  }

  const auto& images = get<0>(m_images);
  const auto& corners = get<1>(m_corners);
  auto& preimages = get<1>(m_images);
  auto& success = get<1>(m_success);

  std::array<Coordinate, dim> v;
  for (unsigned i = 0; i < dim-1; ++i) {
    v[i] = interpolate(images[i+1], corners);
    v[i] -= interpolate(images[0], corners);
  }

  Matrix m;
  for (unsigned i = 0; i < dim-1; ++i) {
    for (unsigned j = 0; j < dim-1; ++j) {
      m[i][j] = v[i]*v[j];
    }
  }

  for (unsigned i = 0; i < dim; ++i) {
    /* Convert yᵢ to barycentric coordinates with respect to Φ(xⱼ) */
    v[dim-1] = corners[i];
    v[dim-1] -= interpolate(images[0], corners);

    Vector rhs, z;
    for (unsigned j = 0; j < dim-1; ++j)
      rhs[j] = v[dim-1]*v[j];
    m.solve(z, rhs);

    //std::cout << "m = " << std::endl << m << std::endl;
    //std::cout << "rhs = " << rhs << std::endl;

    for (unsigned j = 0; j < dim-1; ++j)
      preimages[i][j] = z[j];

    /* Calculate distance along normal direction */
    {
      const auto x = interpolate(z, get<0>(m_corners));
      auto nx = interpolate(z, get<0>(m_normals));
      nx /= nx.two_norm();
      preimages[i][dim-1] = (x - corners[i])*nx;
    }

    /* Check y_i lies inside the Φ(xⱼ) */
    bool success_i = true;
    success_i = success_i && preimages[i][dim-1] <= m_overlap+m_epsilon;
    Field sum(0);
    for (unsigned j = 0; j < dim-1; ++j) {
      success_i = success_i && z[j] >= -m_epsilon;
      sum += z[j];
    }
    success_i = success_i && sum <= 1 + m_epsilon;

    success.set(i, success_i);
  }

#if 0
  std::cout << "-------------------------------------------" << std::endl
            << "inverse projection: success = " << success << std::endl
            << std::endl;
#endif
}

template<typename Coordinate>
void
Projection<Coordinate>
::doEdgeIntersection()
{
  using namespace ProjectionImplementation;
  using std::get;

  m_number_of_edge_intersections = 0;

  /* There are no edge intersections for 2d, only for 3d */
  if (dim != 3)
    return;

  /* There are no edge intersections
   * - when the projection is invalid,
   * - when the projected triangle lies fully in the target triangle,
   * - or when the target triangle lies fully in the projected triangle.
   */
  if (!m_projection_valid || get<0>(m_success).all() || get<1>(m_success).all()) {
    return;
  }

  const auto& images = get<0>(m_images);
  const auto& ys = get<1>(m_corners);

  /* Intersect line through Φ(xᵢ), Φ(xⱼ) with line through yₖ, yₗ:
     We want α, β ∈ ℝ such that
       Φ(xᵢ) + α (Φ(xⱼ) - Φ(xᵢ)) = yₖ + β (yₗ - yₖ)
     or
       α (Φ(xⱼ)-Φ(xᵢ)) + β (yₗ-yₖ) = yₖ-Φ(xᵢ)
     To get a 2×2 system of equations, multiply with yₘ-y₀ for
     m ∈ {1,̣̣2} which are linear indep. (and so the system is
     equivalent to the original 3×2 system)
  */
  for (unsigned edgex = 0; edgex < dim; ++edgex) {
    unsigned i, j;
    std::tie(i, j) = edgeToCorners(edgex);

    /* Both sides of edgex lie in the target triangle means no edge intersection */
    if (m_success.first[i] && m_success.first[j])
      continue;

    const auto pxi = interpolate(images[i], ys);
    const auto pxj = interpolate(images[j], ys);
    const auto pxjpxi = pxj - pxi;

    typedef Dune::FieldMatrix<Field, dim-1, dim-1> Matrix;
    typedef Dune::FieldVector<Field, dim-1> Vector;

    for (unsigned edgey = 0; edgey < dim; ++edgey) {
      unsigned k, l;
      std::tie(k, l) = edgeToCorners(edgey);

      /* Both sides of edgey lie in the projected triangle means no edge intersection */
      if (m_success.second[k] && m_success.second[l])
        continue;

      const auto ykyl = ys[k] - ys[l];
      const auto ykpxi = ys[k] - pxi;

      Matrix mat;
      Vector rhs, z;

      for (unsigned m = 0; m < dim-1; ++m) {
        const auto ym1y0 = ys[m+1] - ys[0];
        mat[m][0] = pxjpxi * ym1y0;
        mat[m][1] = ykyl * ym1y0;
        rhs[m] = ykpxi * ym1y0;
      }

      try {
        mat.solve(z, rhs);

        /* If solving the system gives a NaN, the edges are probably parallel. */
        if (z[0] != z[0] || z[1] != z[1])
          continue;

        Coordinate local_x = corner<Coordinate, Field>(i);
        local_x.axpy(z[0], corner<Coordinate, Field>(j) - corner<Coordinate, Field>(i));
        Coordinate local_y = corner<Coordinate, Field>(k);
        local_y.axpy(z[1], corner<Coordinate, Field>(l) - corner<Coordinate, Field>(k));

#if 0
        std::cout 
                  << "matrix: " << std::endl
                  << mat << std::endl
                  << "rhs: " << rhs << std::endl
                  << "solution: " << z << std::endl
                  << "local_x: " << local_x << std::endl
                  << "local_y: " << local_y << std::endl;

        // TODO: Use earlier check
        if (z[0] != z[0] || z[1] != z[1]) {
          std::cout << "NAN" << std::endl;
          continue;
        }
#endif

        if (!inside(local_x, m_epsilon) || !inside(local_y, m_epsilon))
          continue;

        /* TODO: check overlap */

        auto& intersection = m_edge_intersections[m_number_of_edge_intersections++];
        intersection = { {edgex, edgey}, {local_x, local_y} };

      }
      catch(const Dune::FMatrixError&) {
        /* Edges might be parallel, ignore and continue with next edge */
      }
    }
  }
}

template<typename Coordinate>
bool Projection<Coordinate>
::projectionFeasible(const Coordinate& x, const Coordinate& nx, const Coordinate& y, const Coordinate& ny) const
{
  return true;
}

template<typename Coordinate>
void Projection<Coordinate>
::project()
{
  doProjection();
  doInverseProjection();
  doEdgeIntersection();
}

} /* namespace GridGlue */
} /* namespace Dune */
