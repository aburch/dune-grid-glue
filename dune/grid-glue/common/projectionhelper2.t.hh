#include <dune/common/fmatrix.hh>
#include <limits>

namespace Dune {
namespace GridGlue {

namespace ProjectionImplementation {

template<typename Coordinate, typename Field>
std::pair<bool, unsigned>
isCorner(const Coordinate& x, const Field& epsilon)
{
  using std::abs;
  const static unsigned invalid = std::numeric_limits<unsigned>::max();
  const static unsigned dim = Coordinate::dimension;

  Coordinate tmp = x;
  tmp[dim-1] = Field(0);

  if (tmp.infinity_norm() <= epsilon)
    return {true, 0};

  if (tmp.infinity_norm() < Field(1) - epsilon)
    return {false, invalid};

  for (unsigned i = 0; i < dim - 1; ++i)
    if (abs(x[i]) > Field(1) - epsilon)
      return {true, i + 1};

  return {false, invalid};
}

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

/* Translate edge to corner numbers */
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

template<typename Coordinate>
inline Coordinate
barycentricToCartesian(const Coordinate& x, const std::array<Coordinate, Coordinate::dimension>& corners)
{
  Coordinate y = corners[0];
  for (unsigned i = 0; i < Coordinate::dimension - 1; ++i)
    y.axpy(x[i], corners[i+1] - corners[0]);
  return y;
}

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

} /* namespace ProjectionImplementation */

template<typename Coordinate>
Projection<Coordinate>
::Projection(const std::pair<Corners, Corners>& corners, const std::pair<Normals, Normals>& normals)
  : m_corners(corners)
  , m_normals(normals)
{
  /* Nothing. */
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

  /* Project x_i */
  for (unsigned i = 0; i < origin.size(); ++i) {
    std::cout << "do projection of " << origin[i] << std::endl;
    for (unsigned j = 0; j < dim; ++j)
      m[j][dim-1] = normals[i][j];

    const Coordinate rhs = origin[i] - target[0];
    try {
      std::cout << "m = " << m << std::endl;
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
    v[i] = barycentricToCartesian(images[i+1], corners);
    v[i] -= barycentricToCartesian(images[0], corners);
  }

  Matrix m;
  for (unsigned i = 0; i < dim-1; ++i) {
    for (unsigned j = 0; j < dim-1; ++j) {
      m[i][j] = v[i]*v[j];
    }
  }

  for (unsigned i = 0; i < dim; ++i) {
    /* Convert y_i to barycentric coordinates with respect to \phi(x_j) */
    v[dim-1] = corners[i];
    v[dim-1] -= barycentricToCartesian(images[0], corners);

    Vector rhs, z;
    for (unsigned j = 0; j < dim-1; ++j)
      rhs[j] = v[dim-1]*v[j];
    m.solve(z, rhs);

    std::cout << "m = " << std::endl << m << std::endl;
    std::cout << "rhs = " << rhs << std::endl;

    for (unsigned j = 0; j < dim-1; ++j)
      preimages[i][j] = z[j];

    /* Calculate distance along normal direction */
    {
      const auto& xs = get<0>(m_corners);
      const auto& ns = get<0>(m_normals);
      auto x = xs[0];
      auto nx = ns[0];
      for (unsigned j = 0; j < dim-1; ++j) {
        x.axpy(z[j], xs[j+1] - xs[0]);
        nx.axpy(z[j], ns[j+1] - ns[0]);
      }
      nx /= nx.two_norm();
      preimages[i][dim-1] = (x - corners[i])*nx;
    }

    /* Check y_i lies inside the \phi(x_j) */
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
#if 0
  if (!m_projection_valid || get<0>(m_success).all() || get<1>(m_success).all()) {
    for (auto& in : m_edge_intersections)
      get<0>(in) = false;
    return;
  }
#endif

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
    const auto pxi = barycentricToCartesian(images[i], ys);
    const auto pxj = barycentricToCartesian(images[j], ys);
    const auto pxipxj = pxi - pxj;
    const auto corner_pxi = isCorner(images[i], m_epsilon).second;
    const auto corner_pxj = isCorner(images[j], m_epsilon).second;
    std::cout << "images[i] = " << images[i] << " pxi = " << pxi << " corner_pxi = " << corner_pxi << std::endl;
    std::cout << "images[j] = " << images[j] << " pxj = " << pxj << " corner_pxj = " << corner_pxj << std::endl;

    typedef Dune::FieldMatrix<Field, dim-1, dim-1> Matrix;
    typedef Dune::FieldVector<Field, dim-1> Vector;

    for (unsigned edgey = 0; edgey < dim; ++edgey) {
      unsigned k, l;
      std::tie(k, l) = edgeToCorners(edgey);
      if (corner_pxi == k || corner_pxi == l || corner_pxj == k || corner_pxj == l)
        continue;

      const auto ykyl = ys[k] - ys[l];
      const auto ykpxi = ys[k] - pxi;

      Matrix mat;
      Vector rhs, z;

      for (unsigned m = 0; m < dim-1; ++m) {
        const auto ym1y0 = ys[m+1] - ys[0];
        mat[m][0] = pxipxj * ym1y0;
        mat[m][1] = ykyl * ym1y0;
        rhs[m] = ykpxi * ym1y0;
      }

      try {
        mat.solve(z, rhs);

        Coordinate local_x = corner<Coordinate, Field>(i);
        local_x.axpy(z[0], corner<Coordinate, Field>(j) - corner<Coordinate, Field>(i));
        Coordinate local_y = corner<Coordinate, Field>(k);
        local_y.axpy(z[1], corner<Coordinate, Field>(l) - corner<Coordinate, Field>(k));

        //if (!inside(local_x, m_epsilon) || !inside(local_y, m_epsilon))
        //  continue;

        /* TODO: check overlap */

        std::cout << "-----------------------------" << std::endl
                  << "intersecting lines " << edgex << " and " << edgey << ":" << std::endl
                  << "matrix: " << std::endl
                  << mat << std::endl
                  << "rhs: " << rhs << std::endl
                  << "solution: " << z << std::endl
                  << "local_x: " << local_x << std::endl
                  << "local_y: " << local_y << std::endl;
      }
      catch(const Dune::FMatrixError&) {
        /* Edges might be parallel, ignore and continue with next edge */
      }
    }
  }
}

} /* namespace GridGlue */
} /* namespace Dune */
