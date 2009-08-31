// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SimplicialMeshExtractor.hh
 *  Version:     1.0
 *  Created on:  Feb 19, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for "flat" simplicial grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief
 */

#ifndef SIMPLICIALMESHEXTRACTOR_HH_
#define SIMPLICIALMESHEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/geometry.hh>

#include "surfacedescriptor.hh"
#include <dune/glue/extractors/codim0extractor.hh>


/**
 * @brief grid extractor implementation for simplicial grids
 *
 * Provides methods that build topology information for given grids.
 * Note that these methods only operate on the grid.
 *
 * \tparam GV the grid view class type
 */
template<typename GV>
class SimplicialMeshExtractor
  : public Codim0Extractor<GV>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */


  enum {dim      = GV::dimension};

  enum {dimworld = GV::dimensionworld};

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = dim+1
  };

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                                Coords;
  typedef Dune::array<unsigned int, simplex_corners>                        SimplexTopology;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;
  typedef typename GV::Traits::template Codim<dim>::Iterator VertexIter;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

private:

  /************************** PRIVATE SUBCLASSES **********************/

  /**
   * @class FaceInfo
   * @brief simple struct holding some packed information about this codim 0 entity
   */
  struct FaceInfo
  {
    FaceInfo()
    {}

    FaceInfo(unsigned int index_, IndexType elmtindex_)
      :       self(elmtindex_), index(index_)
    {}

    /// @brief the index of this element (from index set)
    IndexType self;

    /// @brief the index of this face (in internal storage scheme) // NEEDED??
    unsigned int index;

    /// @brief the vertex indices of the corners
    unsigned int corners[simplex_corners];
  };


  /************************** MEMBER VARIABLES ************************/

  // these values are filled on surface extraction and can be
  // asked by the corresponding getters

  /*        Geometrical and Topological Information                */

  /// @brief all information about the extracted faces
  std::vector<FaceInfo>         _faces;


public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  SimplicialMeshExtractor(const GV& gv)
    : Codim0Extractor<GV>(gv)
  {
    std::cout << "This is SimplicialMeshExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid!" << std::endl;
  }


  /**
   * @brief default destructor, frees memory
   */
  ~SimplicialMeshExtractor();


  /*  F U N C T I O N A L I T Y  */

  /**
   */
  void update(const ElementDescriptor<GV>& descr);


  /**
   * @brief delete everything build up so far and free the memory
   */
  void clear();


  /*  S E T T E R S  */


  /*  G E T T E R S  */

  /**
   * @brief getter for the indices array
   * It is strongly recommended not to modify its contents.
   * Deallocation is done in this class.
   * @return the _indices array
   */
  void getFaces(std::vector<SimplexTopology>& faces) const
  {
    faces.resize(this->_faces.size());
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
      for (unsigned int j = 0; j < simplex_corners; ++j)
        faces[i][j] = this->_faces[i].corners[j];
  }


  /**
   * @brief gets index of coordinate in _coords associated with given vertex
   * @return the index if possible, -1 else
   */
  int coordinateIndex(const Vertex& v) const
  {
    typename Codim0Extractor<GV>::VertexInfoMap::const_iterator it = this->_vtxInfo.find(this->index<dim>(v));
    return it == this->_vtxInfo.end() ? -1 : it->second->idx;
  }


  /**
   * @brief gets index of first face as well as the total number of faces that
   * were extracted from this element
   * @param e the element
   * @param first will contain the first index if found, else -1
   * @param count will be 1 if found, else 0.
   * Note: parameter is kept only to have the same interface as the surface extractor
   * @return success
   */
  bool faceIndices(const Element& e, int& first, int& count) const
  {
    typename Codim0Extractor<GV>::ElementInfoMap::const_iterator it = this->_elmtInfo.find(this->indexSet().template index<0>(e));
    if (it == this->_elmtInfo.end())
    {
      first = -1;
      count = 0;
      return false;
    }
    // the iterator is valid, fill the out params
    first = it->second->idx;
    count = 1;
    return true;
  }


  /**
   * @brief gets the number of the face in the element (which is the face)
   *
   * Note: this method is only kept to have the same interface as the surface extractor
   * @param index the index of the element face
   * @return if failed -1, else 0 (consistent with Dune speaking of entity<0>(0) )
   */
  int indexInInside(unsigned int index) const
  {
    return index < this->_faces.size() ? 0 : -1;
  }


  void globalCoords(unsigned int index, const Dune::FieldVector<ctype, dimworld> &bcoords, Coords &wcoords) const;


  void localCoords(unsigned int index, const Dune::FieldVector<ctype, dimworld> &bcoords, Coords &ecoords) const;

  /**
   * @brief for given barycentric coords in a cube compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index, const Dune::FieldVector<ctype, dimworld> &bcoords, Coords &ecoords, Coords &wcoords) const;


  template<typename CoordContainerB, typename CoordContainerW>
  void globalCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerW &wcoords, int size) const;


  template<typename CoordContainerB, typename CoordContainerE>
  void localCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, int size) const;


  /**
   * @brief for several given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index,
                            const Dune::array<Dune::FieldVector<ctype,dim>, dimworld> &subEntityCoords,
                            Dune::array<Dune::FieldVector<ctype,dim>, dimworld> &elementCoords,
                            Dune::array<Dune::FieldVector<ctype,dimworld>, dimworld> &wcoords,
                            int size) const;


  /**
   * @brief gets for each vertex corner of given face (by index) the number of
   * the vertex in parent element's ordering
   * @param index the face's index
   * @param corner the index of the corner
   * @return -1 <=> index invalid or array not filled, else index
   */
  int numCornerInParent(unsigned int index, unsigned int corner) const
  {
    return (index >= this->_faces.size() || corner >= simplex_corners) ?
           -1 : corner;
  }


  /**
   * @brief gets the the element for a given face/element index,
   * throws an exception if index not valid
   * @param index the index of the face
   * @return a reference to the element's stored pointer
   */
  const ElementPtr& element(unsigned int index) const
  {
    if (index >= this->_faces.size())
      DUNE_THROW(Dune::GridError, "invalid face index");
    return (this->_elmtInfo.find(this->_faces[index].self))->second->p;
  }


  /**
   * @brief gets the vertex for a given coordinate index
   * throws an exception if index not valid
   * @param index the index of the coordinate
   * @return a reference to the vertex' stored pointer
   */
  const VertexPtr& vertex(unsigned int index) const
  {
    if (index >= this->_coords.size())
      DUNE_THROW(Dune::GridError, "invalid coordinate index");
    return (this->_vtxInfo.find(this->_coords[index].self))->second->p;
  }


  /**
   * @brief gets the indices of all faces with the given coordinate as corner
   * @param index the index of the coordinate
   * @param faces array if given index was legal
   * @param count length of array if successful
   * @return TRUE <=> if successful
   * DO NOT MODIFY THE ARRAY'S CONTENT!
   */
  bool parentFaces(unsigned int index,  unsigned int const*& faces, unsigned int& count) const
  {
    if (index >= this->_coords.size())
      return false;
    // index valid
    faces = this->_coords[index].faces;
    count = this->_coords[index].num_faces;
    return true;
  }

}; // end of class SimplicialMeshExtractor



template<typename GV>
SimplicialMeshExtractor<GV>::~SimplicialMeshExtractor()
{
  // only the objects that have been allocated manually have to be
  // deallocated manually again
  // free all the manually allocated memory
  for (unsigned int i = 0; i < this->_coords.size(); ++i)
    if (this->_coords[i].faces != NULL)
      delete this->_coords[i].faces;
  for (typename Codim0Extractor<GV>::VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename Codim0Extractor<GV>::ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
}



template<typename GV>
void SimplicialMeshExtractor<GV>::clear()
{
  // this is an inofficial way on how to free the memory allocated
  // by a std::vector
  {
    // free all the manually allocated memory
    for (unsigned int i = 0; i < this->_coords.size(); ++i)
      if (this->_coords[i].faces != NULL)
        delete this->_coords[i].faces;
    std::vector<typename Codim0Extractor<GV>::CoordinateInfo> dummy;
    this->_coords.swap(dummy);
  }
  {
    std::vector<FaceInfo> dummy;
    this->_faces.swap(dummy);
  }

  // first free all manually allocated vertex/element info items...
  for (typename Codim0Extractor<GV>::VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename Codim0Extractor<GV>::ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  // ...then clear the maps themselves, too
  this->_vtxInfo.clear();
  this->_elmtInfo.clear();
}



template<typename GV>
void SimplicialMeshExtractor<GV>::update(const ElementDescriptor<GV>& descr)
{
  // free everything there is in this object
  this->clear();

  // In this first pass iterate over all entities of codim 0.
  // For each codim 1 intersection check if it is part of the boundary and if so,
  // get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the boundary faces' corners.

  {
    // several counter for consecutive indexing are needed
    int simplex_index = 0;
    int vertex_index = 0;
    IndexType eindex = 0;             // supress warning

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<FaceInfo> temp_faces;

    // iterate over all codim 0 elemets on the grid
    for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
    {
      // check if there are unwanted geometric shapes
      // if one appears => exit with error
      if (!elit->type().isSimplex())
        DUNE_THROW(Dune::GridError, "Expected simplicial grid but found a " << elit->type());

      // only do sth. if this element is "interesting"
      // implicit cast is done automatically
      if (descr.contains(elit))
      {
        // add an entry to the element info map, the index will be set properly later
        eindex = this->indexSet().template index<0>(*elit);
        this->_elmtInfo[eindex] = new typename Codim0Extractor<GV>::ElementInfo(simplex_index, elit);


        // add a new face to the temporary collection
        temp_faces.push_back(FaceInfo(simplex_index, eindex));

        // try for each of the faces vertices whether it is already inserted or not
        for (int i = 0; i < simplex_corners; ++i)
        {
          // TODO Think about a way to ensure that orientation of all faces is consistent!
          // get the vertex pointer and the index from the index set
          // Note that the orientation is always the same for all simplices,
          // i.e. CCW which is 0,1 in 2D and 0,1,2 in 3D
          typename Codim0Extractor<GV>::VertexPtr vptr(elit->template subEntity<dim>(i));
          IndexType vindex = this->indexSet().template index<dim>(*vptr);

          // if the vertex is not yet inserted in the vertex info map
          // it is a new one -> it will be inserted now!
          typename Codim0Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
          if (vimit == this->_vtxInfo.end())
          {
            // insert into the map
            this->_vtxInfo[vindex] = new typename Codim0Extractor<GV>::VertexInfo(vertex_index, vptr);
            // remember the vertex as a corner of the current face in temp_faces
            temp_faces.back().corners[i] = vertex_index;
            // increase the current index
            vertex_index++;
          }
          else
          {
            // only insert the index into the simplices array
            temp_faces.back().corners[i] = vimit->second->idx;
          }
        }

        // now increase the current face index
        simplex_index++;
      }
    }             // end loop over elements

    // allocate the array for the face specific information...
    this->_faces.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());
  }


  // now first write the array with the coordinates...
  this->_coords.resize(this->_vtxInfo.size());
  typename Codim0Extractor<GV>::VertexInfoMap::const_iterator it1 = this->_vtxInfo.begin();
  for (; it1 != this->_vtxInfo.end(); ++it1)
  {
    // get a pointer to the associated info object
    typename Codim0Extractor<GV>::CoordinateInfo* current = &this->_coords[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->self = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }


  // now add the vertices' parent faces in the _vertexFaces map.
  // therefore iterate over all indices in the _index array...
  {
    std::vector<unsigned int> refcount(this->_coords.size(), 0);

    // for each coordinate count the references in the _indices array
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
      for (unsigned int j = 0; j < simplex_corners; ++j)
        refcount[this->_faces[i].corners[j]]++;

    // allocate the right amount of storage for each vertex's references
    for (unsigned int i = 0; i < this->_coords.size(); ++i)
    {
      // allocate an array and initialize its first element with its length
      this->_coords[i].num_faces = refcount[i];
      this->_coords[i].faces = new unsigned int[refcount[i]];
      refcount[i] = 0;                   // used as "pointer" in next loop
    }

    // add the references
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
    {
      for (unsigned int j = 0; j < simplex_corners; ++j)
      {
        unsigned int ref = this->_faces[i].corners[j];
        this->_coords[ref].faces[refcount[ref]] = i;
        refcount[ref]++;
      }
    }
  }


  //	const char prefix[] = "SimplicialMeshExtractor: ";
  //
  //	STDOUTLN(prefix << "Extracted Coordinates (size=" << this->_coords.size() << ")");
  //	for (unsigned int i = 0; i < this->_coords.size(); ++i)
  //	{
  ////		if (i % 100 == 0)
  //		{
  //			STDOUT(prefix << "self=" << this->_coords[i].self << " index=" << this->_coords[i].index
  //					<< " coord=(" << this->_coords[i].coord << ") num_faces=" << this->_coords[i].num_faces
  //					<< " faces={");
  //			for (unsigned int j = 0; j < this->_coords[i].num_faces; ++j)
  //				STDOUT(" " << this->_coords[i].faces[j]);
  //			STDOUTLN("}");
  //		}
  //	}
  //
  //	STDOUTLN("\n" << prefix << "Extracted faces (size=" << this->_faces.size() << ")");
  //	for (unsigned int i = 0; i < this->_faces.size(); ++i)
  //	{
  ////		if (i % 100 == 0)
  //		{
  //			STDOUT(prefix << "self=" << this->_faces[i].self << " index=" << this->_faces[i].index << " corners={");
  //			for (unsigned int j = 0; j < simplex_corners; ++j)
  //				STDOUT("(" << this->_faces[i].corners[j] << ")");
  //			STDOUTLN("}");
  //		}
  //
  //	}

}


template<typename GV>
inline void SimplicialMeshExtractor<GV>::globalCoords(unsigned int index, const Dune::FieldVector<ctype, dimworld> &bcoords, Coords &wcoords) const
{
  wcoords = this->_elmtInfo.find(this->_faces[index].self)->second->p->geometry().global(barycentricToReference(bcoords));
}


template<typename GV>
inline void SimplicialMeshExtractor<GV>::localCoords(unsigned int index, const Dune::FieldVector<ctype, dimworld> &bcoords, Coords &ecoords) const
{
  ecoords = barycentricToReference(bcoords);
}


template<typename GV>
inline void SimplicialMeshExtractor<GV>::localAndGlobalCoords(unsigned int index, const Dune::FieldVector<ctype, dimworld> &bcoords, Coords &ecoords, Coords &wcoords) const
{
  this->localCoords(index, bcoords, ecoords);
  this->globalCoords(index, bcoords, wcoords);
}


template<typename GV>
template<typename CoordContainerB, typename CoordContainerW>
void SimplicialMeshExtractor<GV>::globalCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerW &wcoords, int size) const
{
  ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].self)->second->p;
  for (int i = 0; i < size; ++i)
  {
    // compute global coordinates and pad with 0 in last coordinate
    wcoords[i] = eptr->geometry().global(barycentricToReference(bcoords[i]));
  }
}


template<typename GV>
template<typename CoordContainerB, typename CoordContainerE>
void SimplicialMeshExtractor<GV>::localCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, int size) const
{
  for (int i = 0; i < size; ++i)
    ecoords[i] = barycentricToReference(bcoords[i]);
}


template<typename GV>
void SimplicialMeshExtractor<GV>::
localAndGlobalCoords(unsigned int index,
                     const Dune::array<Dune::FieldVector<ctype,dim>, dimworld> &subEntityCoords,
                     Dune::array<Dune::FieldVector<ctype,dim>, dimworld> &elementCoords,
                     Dune::array<Dune::FieldVector<ctype,dimworld>, dimworld> &wcoords,
                     int size) const
{
  ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].self)->second->p;
  for (int i = 0; i < size; ++i)
  {
    elementCoords[i] = subEntityCoords[i];
    // compute global coordinates
    wcoords[i] = eptr->geometry().global(elementCoords[i]);
  }
}



// TODO add all those local and global coordinates functions

#endif // SIMPLICIALMESHEXTRACTOR_HH_