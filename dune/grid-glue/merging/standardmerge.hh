// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * \brief Common base class for many merger implementations: produce pairs of entities that _may_ intersect
 */

#ifndef STANDARD_MERGE_HH
#define STANDARD_MERGE_HH


#include <iostream>
#include <iomanip>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/version.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/merger.hh>



/** \brief Common base class for many merger implementations: produce pairs of entities that _may_ intersect

   Many merger algorithms consist of two parts: on the one hand there is a mechanism that produces pairs of
   elements that may intersect.  On the other hand there is an algorithm that computes the intersection of two
   given elements.  For the pairs-producing algorithm there appears to be a canonical choice, namely the algorithm
   by Gander and Japhet described in 'An Algorithm for Non-Matching Grid Projections with Linear Complexity,
   M.J. Gander and C. Japhet, Domain Decomposition Methods in Science and Engineering XVIII,
   pp. 185--192, Springer-Verlag, 2009.'  This class implements this algorithm, calling a pure virtual function
   computeIntersection() to compute the intersection between two elements.  Actual merger implementations
   can derive from this class and only implement computeIntersection().

   \tparam T The type used for coordinates (assumed to be the same for both grids)
   \tparam grid1Dim Dimension of the grid1 grid
   \tparam grid2Dim Dimension of the grid2 grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template<class T, int grid1Dim, int grid2Dim, int dimworld>
class StandardMerge
  : public Merger<T,grid1Dim,grid2Dim,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief Type used for local coordinates on the grid1 side
  typedef typename Merger<T,grid1Dim,grid2Dim,dimworld>::Grid1Coords Grid1Coords;

  /// @brief Type used for local coordinates on the grid2 side
  typedef typename Merger<T,grid1Dim,grid2Dim,dimworld>::Grid2Coords Grid2Coords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

protected:

  bool valid;

  StandardMerge() : valid(false) {}

  struct RemoteSimplicialIntersection
  {
    /** \brief Dimension of this intersection */
    enum {intersectionDim = (grid1Dim<grid2Dim) ? grid1Dim : grid2Dim};

    /** \brief Number of vertices of the intersection (it's a simplex) */
    enum {nVertices = intersectionDim + 1};

    /** \brief Default constructor */
    RemoteSimplicialIntersection()
    {}

    /** \brief Constructor for two given entity indices */
    RemoteSimplicialIntersection(int grid1Entity, int grid2Entity)
    : grid1Entity_(grid1Entity), grid2Entity_(grid2Entity)
    {}

    // Local coordinates in the grid1 entity
    Dune::array<Dune::FieldVector<T,grid1Dim>, nVertices> grid1Local_;

    // Local coordinates in the grid1 entity
    Dune::array<Dune::FieldVector<T,grid2Dim>, nVertices> grid2Local_;

    //
    int grid1Entity_;

    int grid2Entity_;

  };

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  virtual void computeIntersection(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   unsigned int grid1Index,
                                   std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   unsigned int grid2Index,
                                   std::bitset<(1<<grid2Dim)>& neighborIntersects2) = 0;

  /** \brief Compute the intersection between two overlapping elements
   * \return true if there was a nonempty intersection
   */
  bool computeIntersection(unsigned int candidate0, unsigned int candidate1,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                           const std::vector<Dune::GeometryType>& grid1_element_types,
                           std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                           const std::vector<Dune::GeometryType>& grid2_element_types,
                           std::bitset<(1<<grid2Dim)>& neighborIntersects2);

  /** \brief Test whether there is a nonempty intersection between two overlapping elements
   * \return true if there was a nonempty intersection
   */
  bool testIntersection(unsigned int candidate0, unsigned int candidate1,
                        const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                        const std::vector<Dune::GeometryType>& grid1_element_types,
                        std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                        const std::vector<Dune::GeometryType>& grid2_element_types,
                        std::bitset<(1<<grid2Dim)>& neighborIntersects2);

  int bruteForceSearch(int candidate1,
                       const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                       const std::vector<Dune::GeometryType>& grid1_element_types,
                       const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                       const std::vector<Dune::GeometryType>& grid2_element_types);

  template <int gridDim>
  void computeNeighborsPerElement(const std::vector<Dune::GeometryType>& gridElementTypes,
                                  const std::vector<std::vector<unsigned int> >& gridElementCorners,
                                  std::vector<std::vector<int> >& elementNeighbors);

  /*   M E M B E R   V A R I A B L E S   */

  /** \brief The computed intersections */
  std::vector<RemoteSimplicialIntersection> intersections_;

  /** \brief Temporary internal data */
  std::vector<std::vector<unsigned int> > grid1ElementCorners_;
  std::vector<std::vector<unsigned int> > grid2ElementCorners_;

  std::vector<std::vector<int> > elementNeighbors1_;
  std::vector<std::vector<int> > elementNeighbors2_;

public:

  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @copydoc Merger<T,grid1Dim,grid2Dim,dimworld>::build
   */
  virtual void build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1_Coords,
             const std::vector<unsigned int>& grid1_elements,
             const std::vector<Dune::GeometryType>& grid1_element_types,
             const std::vector<Dune::FieldVector<T,dimworld> >& grid2_coords,
             const std::vector<unsigned int>& grid2_elements,
             const std::vector<Dune::GeometryType>& grid2_element_types
             );


  /*   Q U E S T I O N I N G   T H E   M E R G E D   G R I D   */

  /// @brief get the number of simplices in the merged grid
  /// The indices are then in 0..nSimplices()-1
  unsigned int nSimplices() const;

  void clear()
  {
    // Delete old internal data, from a possible previous run
    purge(intersections_);
    purge(grid1ElementCorners_);
    purge(grid2ElementCorners_);

    valid = false;
  }

private:

  /** clear arbitrary containers */
  template<typename V>
  static void purge(V & v)
  {
    v.clear();
    V v2(v);
    v.swap(v2);
  }

  /*   M A P P I N G   O N   I N D E X   B A S I S   */

  /**
   * @brief get index of grid1 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid1 parent simplex
   */
  unsigned int grid1Parent(unsigned int idx) const;

  /**
   * @brief get index of grid2 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid2 parent simplex
   */
  unsigned int grid2Parent(unsigned int idx) const;


  /*   G E O M E T R I C A L   I N F O R M A T I O N   */

  /**
   * @brief get the grid1 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid1Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid1 simplex
   */
  Grid1Coords grid1ParentLocal(unsigned int idx, unsigned int corner) const;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2 simplex
   */
  Grid2Coords grid2ParentLocal(unsigned int idx, unsigned int corner) const;

};


/* IMPLEMENTATION */

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::computeIntersection(unsigned int candidate0, unsigned int candidate1,
                                                                      const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                      const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                      std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                                                      const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                      const std::vector<Dune::GeometryType>& grid2_element_types,
                                                                      std::bitset<(1<<grid2Dim)>& neighborIntersects2)
{
  unsigned int oldNumberOfIntersections = intersections_.size();

  // Select vertices of the grid1 element
  int grid1NumVertices = grid1ElementCorners_[candidate0].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid1ElementCorners(grid1NumVertices);
  for (int i=0; i<grid1NumVertices; i++)
    grid1ElementCorners[i] = grid1Coords[grid1ElementCorners_[candidate0][i]];

  // Select vertices of the grid2 element
  int grid2NumVertices = grid2ElementCorners_[candidate1].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid2ElementCorners(grid2NumVertices);
  for (int i=0; i<grid2NumVertices; i++)
    grid2ElementCorners[i] = grid2Coords[grid2ElementCorners_[candidate1][i]];

  // ///////////////////////////////////////////////////////
  //   Compute the intersection between the two elements
  // ///////////////////////////////////////////////////////

  computeIntersection(grid1_element_types[candidate0], grid1ElementCorners, candidate0, neighborIntersects1,
                      grid2_element_types[candidate1], grid2ElementCorners, candidate1, neighborIntersects2);

  // Have we found an intersection?
  return (intersections_.size() > oldNumberOfIntersections);

}


/** \note The implementation of this is very ugly: we compute the entire intersection,
 * which gets appended to the list of all intersections.  Then we throw it away again.
 * I'm afraid in the interest of speed this will have to be another pure virtual method.
 * The implementation (child) classes can usually determine whether an intersection
 * is nonempty without having to compute it entirely first.
 */
template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::testIntersection(unsigned int candidate0, unsigned int candidate1,
                                                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                   const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                   std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                   const std::vector<Dune::GeometryType>& grid2_element_types,
                                                                   std::bitset<(1<<grid2Dim)>& neighborIntersects2)
{
  unsigned int oldNumberOfIntersections = intersections_.size();

  // Select vertices of the grid1 element
  int grid1NumVertices = grid1ElementCorners_[candidate0].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid1ElementCorners(grid1NumVertices);
  for (int i=0; i<grid1NumVertices; i++)
    grid1ElementCorners[i] = grid1Coords[grid1ElementCorners_[candidate0][i]];

  // Select vertices of the grid2 element
  int grid2NumVertices = grid2ElementCorners_[candidate1].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid2ElementCorners(grid2NumVertices);
  for (int i=0; i<grid2NumVertices; i++)
    grid2ElementCorners[i] = grid2Coords[grid2ElementCorners_[candidate1][i]];

  // ///////////////////////////////////////////////////////
  //   Compute the intersection between the two elements
  // ///////////////////////////////////////////////////////

  computeIntersection(grid1_element_types[candidate0], grid1ElementCorners, candidate0, neighborIntersects1,
                      grid2_element_types[candidate1], grid2ElementCorners, candidate1, neighborIntersects2);

  // Have we found an intersection?
  bool intersectionFound = (intersections_.size() > oldNumberOfIntersections);

  // very stupid: delete the compute intersection again
  while (intersections_.size()> oldNumberOfIntersections)
    intersections_.pop_back();

  return intersectionFound;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::bruteForceSearch(int candidate1,
                                                                  const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                  const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                  const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                  const std::vector<Dune::GeometryType>& grid2_element_types)
{
  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    std::bitset<(1<<grid1Dim)> neighborIntersects1;
    std::bitset<(1<<grid2Dim)> neighborIntersects2;
    bool intersectionFound = testIntersection(i, candidate1,
                                              grid1Coords,grid1_element_types, neighborIntersects1,
                                              grid2Coords,grid2_element_types, neighborIntersects2);

    // if there is an intersection, i is our new seed candidate on the grid1 side
    if (intersectionFound)
      return i;

  }

  return -1;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
template<int gridDim>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::
computeNeighborsPerElement(const std::vector<Dune::GeometryType>& gridElementTypes,
                           const std::vector<std::vector<unsigned int> >& gridElementCorners,
                           std::vector<std::vector<int> >& elementNeighbors)
{
  typedef std::vector<unsigned int> FaceType;
  typedef std::map<FaceType, std::pair<unsigned int, unsigned int> > FaceSetType;

  ///////////////////////////////////////////////////////////////////////////////////////
  //  First: grid 1
  ///////////////////////////////////////////////////////////////////////////////////////
  FaceSetType faces;
  elementNeighbors.resize(gridElementTypes.size());

  for (size_t i=0; i<gridElementTypes.size(); i++)
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    elementNeighbors[i].resize(Dune::ReferenceElements<T,gridDim>::general(gridElementTypes[i]).size(1), -1);
#else
    elementNeighbors[i].resize(Dune::GenericReferenceElements<T,gridDim>::general(gridElementTypes[i]).size(1), -1);
#endif

  for (size_t i=0; i<gridElementTypes.size(); i++) {

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    const Dune::ReferenceElement<T,gridDim>& refElement = Dune::ReferenceElements<T,gridDim>::general(gridElementTypes[i]);
#else
    const Dune::GenericReferenceElement<T,gridDim>& refElement = Dune::GenericReferenceElements<T,gridDim>::general(gridElementTypes[i]);
#endif

    for (size_t j=0; j<(size_t)refElement.size(1); j++) {

      FaceType face;
      // extract element face
      for (size_t k=0; k<(size_t)refElement.size(j,1,gridDim); k++)
        face.push_back(gridElementCorners[i][refElement.subEntity(j,1,k,gridDim)]);

      // sort the face vertices to get rid of twists and other permutations
      std::sort(face.begin(), face.end());

      typename FaceSetType::iterator faceHandle = faces.find(face);

      if (faceHandle == faces.end()) {

        // face has not been visited before
        faces.insert(std::make_pair(face, std::make_pair(i,j)));

      } else {

        // face has been visited before: store the mutual neighbor information
        elementNeighbors[i][j] = faceHandle->second.first;
        elementNeighbors[faceHandle->second.first][faceHandle->second.second] = i;

        faces.erase(faceHandle);

      }

    }

  }
}

// /////////////////////////////////////////////////////////////////////
//   Compute the intersection of all pairs of elements
//   Linear algorithm by Gander and Japhet, Proc. of DD18
// /////////////////////////////////////////////////////////////////////

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                        const std::vector<unsigned int>& grid1_elements,
                                                        const std::vector<Dune::GeometryType>& grid1_element_types,
                                                        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                        const std::vector<unsigned int>& grid2_elements,
                                                        const std::vector<Dune::GeometryType>& grid2_element_types
                                                        )
{

  std::cout << "StandardMerge building merged grid..." << std::endl;
  Dune::Timer watch;

  clear();
  // clear global intersection list
  intersections_.clear();
  this->counter = 0;

  // /////////////////////////////////////////////////////////////////////
  //   Copy element corners into a data structure with block-structure.
  //   This is not as efficient but a lot easier to use.
  //   We may think about efficiency later.
  // /////////////////////////////////////////////////////////////////////

  // first the grid1 side
  grid1ElementCorners_.resize(grid1_element_types.size());

  unsigned int grid1CornerCounter = 0;

  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    // Select vertices of the grid1 element
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    int numVertices = Dune::ReferenceElements<T,grid1Dim>::general(grid1_element_types[i]).size(grid1Dim);
#else
    int numVertices = Dune::GenericReferenceElements<T,grid1Dim>::general(grid1_element_types[i]).size(grid1Dim);
#endif
    grid1ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid1ElementCorners_[i][j] = grid1_elements[grid1CornerCounter++];

  }

  // then the grid2 side
  grid2ElementCorners_.resize(grid2_element_types.size());

  unsigned int grid2CornerCounter = 0;

  for (std::size_t i=0; i<grid2_element_types.size(); i++) {

    // Select vertices of the grid2 element
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    int numVertices = Dune::ReferenceElements<T,grid2Dim>::general(grid2_element_types[i]).size(grid2Dim);
#else
    int numVertices = Dune::GenericReferenceElements<T,grid2Dim>::general(grid2_element_types[i]).size(grid2Dim);
#endif
    grid2ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid2ElementCorners_[i][j] = grid2_elements[grid2CornerCounter++];

  }

  ////////////////////////////////////////////////////////////////////////
  //  Compute the face neighbors for each element
  ////////////////////////////////////////////////////////////////////////

  computeNeighborsPerElement<grid1Dim>(grid1_element_types, grid1ElementCorners_, elementNeighbors1_);
  computeNeighborsPerElement<grid2Dim>(grid2_element_types, grid2ElementCorners_, elementNeighbors2_);

  std::cout << "setup took " << watch.elapsed() << " seconds." << std::endl;

  ////////////////////////////////////////////////////////////////////////
  //   Data structures for the advancing-front algorithm
  ////////////////////////////////////////////////////////////////////////

  std::stack<unsigned int> candidates1;
  std::stack<unsigned int> candidates2;

  std::vector<int> seeds(grid2_element_types.size(), -1);

  // /////////////////////////////////////////////////////////////////////
  //   Do a brute-force search to find one pair of intersecting elements
  //   to start the advancing-front type algorithm with.
  // /////////////////////////////////////////////////////////////////////

  // Set flag if element has been handled
  Dune::BitSetVector<1> isHandled2(grid2_element_types.size());

  // Set flag if the element has been entered in the queue
  Dune::BitSetVector<1> isCandidate2(grid2_element_types.size());

  for (std::size_t j=0; j<grid2_element_types.size(); j++) {

    int seed = bruteForceSearch(j,grid1Coords,grid1_element_types,grid2Coords,grid2_element_types);

    if (seed >=0) {
      candidates2.push(j);        // the candidate and a seed for the candidate
      seeds[j] = seed;
      break;
    } else // If the brute force search did not find any intersection we can skip this element
        isHandled2[j] = true;

  }

  // /////////////////////////////////////////////////////////////////////
  //   Main loop
  // /////////////////////////////////////////////////////////////////////

  std::set<unsigned int> isHandled1;

  std::set<unsigned int> isCandidate1;

  while (!candidates2.empty()) {

    // Get the next element on the grid2 side
    unsigned int currentCandidate2 = candidates2.top();
    unsigned int seed = seeds[currentCandidate2];
    assert(seed >= 0);

    candidates2.pop();
    isHandled2[currentCandidate2] = true;

    // Start advancing front algorithm on the grid1 side from the 'seed' element that
    // we stored along with the current grid2 element
    candidates1.push(seed);

    isHandled1.clear();
    isCandidate1.clear();

    while (!candidates1.empty()) {

      unsigned int currentCandidate1 = candidates1.top();
      candidates1.pop();
      isHandled1.insert(currentCandidate1);

      // Test whether there is an intersection between currentCandidate0 and currentCandidate1
      std::bitset<(1<<grid1Dim)> neighborIntersects1;
      std::bitset<(1<<grid2Dim)> neighborIntersects2;
      bool intersectionFound = computeIntersection(currentCandidate1, currentCandidate2,
                                                   grid1Coords,grid1_element_types, neighborIntersects1,
                                                   grid2Coords,grid2_element_types, neighborIntersects2);

      for (size_t i=0; i<neighborIntersects2.size(); i++)
        if (neighborIntersects2[i] and elementNeighbors2_[currentCandidate2][i] != -1)
          seeds[elementNeighbors2_[currentCandidate2][i]] = currentCandidate1;

      // add neighbors of candidate0 to the list of elements to be checked
      if (intersectionFound) {

        for (size_t i=0; i<elementNeighbors1_[currentCandidate1].size(); i++) {

          int neighbor = elementNeighbors1_[currentCandidate1][i];

          if (neighbor == -1)            // do nothing at the grid boundary
            continue;

          if (isHandled1.find(neighbor) == isHandled1.end()
              && isCandidate1.find(neighbor) == isCandidate1.end()) {
            candidates1.push(neighbor);
            isCandidate1.insert(neighbor);
          }

        }

      }

    }

    // We have now found all intersections of elements in the grid1 side with currentCandidate2
    // Now we add all neighbors of currentCandidate2 that have not been treated yet as new
    // candidates.

    // Do we have an unhandled neighbor with a seed?
    bool seedFound = false;
    for (size_t i=0; i<elementNeighbors2_[currentCandidate2].size(); i++) {

      int neighbor = elementNeighbors2_[currentCandidate2][i];

      if (neighbor == -1)         // do nothing at the grid boundary
        continue;

      // Add all unhandled intersecting neighbors to the queue
      if (!isHandled2[neighbor][0] && !isCandidate2[neighbor][0] and seeds[neighbor]>-1) {

        isCandidate2[neighbor][0] = true;
        candidates2.push(neighbor);
        seedFound = true;
      }

    }

    if (seedFound)
      continue;

    // There is no neighbor with a seed, so we need to be a bit more aggressive...
    // get all neighbors of currentCandidate2, but not currentCandidate2 itself
    for (size_t i=0; i<elementNeighbors2_[currentCandidate2].size(); i++) {

      int neighbor = elementNeighbors2_[currentCandidate2][i];

      if (neighbor == -1)         // do nothing at the grid boundary
        continue;

      if (!isHandled2[neighbor][0] && !isCandidate2[neighbor][0]) {

        // Get a seed element for the new grid2 element
        // Look for an element on the grid1 side that intersects the new grid2 element.
        int seed = -1;

        // Look among the ones that have been tested during the last iteration.
        for (typename std::set<unsigned int>::iterator seedIt = isHandled1.begin();
             seedIt != isHandled1.end(); ++seedIt) {

          std::bitset<(1<<grid1Dim)> neighborIntersects1;
          std::bitset<(1<<grid2Dim)> neighborIntersects2;
          bool intersectionFound = testIntersection(*seedIt, neighbor,
                                                    grid1Coords,grid1_element_types, neighborIntersects1,
                                                    grid2Coords,grid2_element_types, neighborIntersects2);

          // if the intersection is nonempty, *seedIt is our new seed candidate on the grid1 side
          if (intersectionFound) {
            seed = *seedIt;
            break;
          }

        }

        if (seed < 0) {
          // The fast method didn't find a grid1 element that intersects with
          // the new grid2 candidate.  We have to do a brute-force search.
          seed = bruteForceSearch(neighbor,
                                  grid1Coords,grid1_element_types,
                                  grid2Coords,grid2_element_types);

        }

        // We have tried all we could: the candidate is 'handled' now
        isCandidate2[neighbor] = true;

        if (seed < 0)
          // still no seed?  Then the new grid2 candidate isn't overlapped by anything
          continue;

        // we have a seed now
        candidates2.push(neighbor);
        seeds[neighbor] = seed;

      }

    }

  }

  valid = true;
  std::cout << "intersection construction took " << watch.elapsed() << " seconds." << std::endl;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::nSimplices() const
{
  assert(valid);
  return intersections_.size();
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1Parent(unsigned int idx) const
{
  assert(valid);
  return intersections_[idx].grid1Entity_;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2Parent(unsigned int idx) const
{
  assert(valid);
  return intersections_[idx].grid2Entity_;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
typename StandardMerge<T,grid1Dim,grid2Dim,dimworld>::Grid1Coords StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1ParentLocal(unsigned int idx, unsigned int corner) const
{
  assert(valid);
  return intersections_[idx].grid1Local_[corner];
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
typename StandardMerge<T,grid1Dim,grid2Dim,dimworld>::Grid2Coords StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2ParentLocal(unsigned int idx, unsigned int corner) const
{
  assert(valid);
  return intersections_[idx].grid2Local_[corner];
}

#define DECL extern
#define STANDARD_MERGE_INSTANTIATE(T,A,B,C) \
  DECL template \
  void StandardMerge<T,A,B,C>::build(const std::vector<Dune::FieldVector<T,C> >& grid1Coords, \
                                     const std::vector<unsigned int>& grid1_elements, \
                                     const std::vector<Dune::GeometryType>& grid1_element_types, \
                                     const std::vector<Dune::FieldVector<T,C> >& grid2Coords, \
                                     const std::vector<unsigned int>& grid2_elements, \
                                     const std::vector<Dune::GeometryType>& grid2_element_types \
                                     )

STANDARD_MERGE_INSTANTIATE(double,1,1,1);
STANDARD_MERGE_INSTANTIATE(double,2,2,2);
STANDARD_MERGE_INSTANTIATE(double,3,3,3);
#undef STANDARD_MERGE_INSTANTIATE
#undef DECL

#endif // STANDARD_MERGE_HH
