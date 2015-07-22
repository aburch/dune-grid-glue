#include "config.h"

#include <dune/grid-glue/common/projectionhelper2.hh>

using Dune::GridGlue::Projection;

bool
test_project_simple()
{
  bool success = true;

  typedef Dune::FieldVector<double, 3> V;
  typedef Projection<V> P;

  const std::pair<P::Corners, P::Corners> corners = {
    {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
    {{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}}
  };
  const std::pair<P::Normals, P::Normals> normals = {
    {{{0, 0, 1}, {0, 0, 1}, {0, 0, 1}}},
    {{{0, 0, 1}, {0, 0, 1}, {0, 0, 1}}}
  };

  P p(corners, normals);
  p.doProjection();
  p.doInverseProjection();

  std::cout << "projection:" << std::endl;
  for (const auto& px : p.m_images.first) {
    std::cout << px << std::endl;
  }
  std::cout << p.m_success.first << std::endl;
  std::cout << "inverse projection:" << std::endl;
  for (const auto& py : p.m_images.second) {
    std::cout << py << std::endl;
  }
  std::cout << p.m_success.second << std::endl;

  p.doEdgeIntersection();

  return success;
}

bool
test_project_simple2()
{
  bool success = true;

  typedef Dune::FieldVector<double, 3> V;
  typedef Projection<V> P;

  const std::pair<P::Corners, P::Corners> corners = {
    {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
    {{{0, 0, 1}, {2, 0, 1}, {0, 2, 1}}}
  };
  const std::pair<P::Normals, P::Normals> normals = {
    {{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}},
    {{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}}
  };

  P p(corners, normals);
  p.doProjection();
  p.doInverseProjection();

  std::cout << "projection:" << std::endl;
  for (const auto& px : p.m_images.first) {
    std::cout << px << std::endl;
  }
  std::cout << p.m_success.first << std::endl;
  std::cout << "inverse projection:" << std::endl;
  for (const auto& py : p.m_images.second) {
    std::cout << py << std::endl;
  }
  std::cout << p.m_success.second << std::endl;

  p.doEdgeIntersection();

  return success;
}

int main()
{
  bool pass = true;

  pass = pass && test_project_simple();
  pass = pass && test_project_simple2();

  return pass ? 0 : 1;
}
