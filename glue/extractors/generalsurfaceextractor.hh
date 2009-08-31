// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GeneralSurfaceExtractor.hh
 *  Version:     1.0
 *  Created on:  Feb 3, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for hybrid surface grids (tris and/or quads on surface)
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief grid extractor implementation for hybrid surface grids (tris and/or quads on surface)
 */

#ifndef GENERALSURFACEEXTRACTOR_HH_
#define GENERALSURFACEEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/geometry.hh>

#include "surfacedescriptor.hh"
#include <dune/glue/extractors/codim1extractor.hh>

/**
 * @brief grid extractor implementation for arbitrary surface grids (tris and/or quads on surface)
 *
 * Provides methods that build topology information for given grids.
 * \tparam GV the grid class type
 */
template<typename GV>
class GeneralSurfaceExtractor
  : public Codim1Extractor<GV>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum
  {
    dimworld = GV::dimensionworld
  };

  enum
  {
    dim = GV::dimension
  };

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = dim
  };

  enum
  {
    cube_corners = 1 << (dim-1)
  };

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                           Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim1Extractor<GV>::FaceInfo FaceInfo;

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  GeneralSurfaceExtractor(const GV& gv)
    : Codim1Extractor<GV>(gv)
  {
    std::cout << "This is GeneralSurfaceExtractor on a <" << GV::dimension
              << "," << GV::dimensionworld << "> grid working in " << dimworld << " space!" << std::endl;
  }

  /*  F U N C T I O N A L I T Y  */


  /**
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a mapper or set object
   * @c m specifying an index set with codimension 0 entities near or on the boundary.
   * It is assumed that only one geometric shape exists on the boundary.
   * The template parameter n then denotes the number of corners per boundary element
   * (e.g. n==3 for triangles in a 3D grid of tetrahedra).
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are cubes.
   * @param descr a descriptor that "selects" the faces to add to the surface
   */
  void update(const FaceDescriptor<GV>& descr);


  /**
   * @brief for given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  void globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const;


  /**
   * @brief for given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  void localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const;


  /**
   * @brief for given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const;


  /**
   * @brief for several given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  template<typename CoordContainer>
  void globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords) const;


  /**
   * @brief for several given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  template<typename CoordContainer>
  void localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords) const;


  /**
   * @brief for several given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   * @return
   */
  void localAndGlobalCoords(unsigned int index,
                            const Dune::array<Dune::FieldVector<ctype,dim-1>, dimworld> &bcoords,
                            Dune::array<Dune::FieldVector<ctype,dim>, dimworld> &ecoords,
                            Dune::array<Dune::FieldVector<ctype,dimworld>, dimworld> &wcoords) const;

}; // end of class GeneralSurfaceExtractor



template<typename GV>
void GeneralSurfaceExtractor<GV>::update(const FaceDescriptor<GV>& descr)
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

    // needed later for insertion into a std::set which only
    // works with const references

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<FaceInfo> temp_faces;

    // iterate over all codim 0 elemets on the grid
    for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
    {
      Dune::GeometryType gt = elit->type();

      // remember the indices of the faces that shall become
      // part of the surface
      std::set<int> boundary_faces;

      // iterate over all intersections of codim 1 and test if the
      // boundary intersections are to be added to the surface
      for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
      {
        // only look at boundary faces
        if (is->boundary() && descr.contains(elit, is->indexInInside()))
          boundary_faces.insert(is->indexInInside());
      }

      // if some face is part of the surface add it!
      if (boundary_faces.size() != 0)
      {
        // add an entry to the element info map, the index will be set properly later,
        // whereas the number of faces is already known
        eindex = this->indexSet().template index<0>(*elit);
        this->_elmtInfo[eindex] = new typename Codim1Extractor<GV>::ElementInfo(simplex_index, elit, 0);

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (typename std::set<int>::const_iterator sit = boundary_faces.begin(); sit != boundary_faces.end(); ++sit)
        {
          // get the corner count of this face
          const int face_corners = Dune::GenericReferenceElements<ctype, dim>::general(gt).size(*sit, 1, dim);

          // register the additional face(s)
          this->_elmtInfo[eindex]->faces += (face_corners - 2);

          // now we only have to care about the 3D case, i.e. a triangle face can be
          // inserted directly whereas a quadrilateral face has to be divided into two triangles
          switch (face_corners)
          {
          case 2 :
            assert(dim == 2);
            // we have a triangle here

            // add a new face to the temporary collection
            temp_faces.push_back(FaceInfo(simplex_index, eindex, *sit, FaceInfo::triangle));

            // try for each of the faces vertices whether it is already inserted or not
            for (int i = 0; i < face_corners; ++i)
            {
              // get the number of the vertex in the parent element
              int vertex_number = orientedSubface<2>(gt, *sit, i);

              // get the vertex pointer and the index from the index set
              VertexPtr vptr(elit->template subEntity<dim>(vertex_number));
              IndexType vindex = this->indexSet().template index<dim>(*vptr);

              // remember the vertex' number in parent element's vertices
              temp_faces.back().corners[i].num = vertex_number;

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename Codim1Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
              if (vimit == this->_vtxInfo.end())
              {
                // insert into the map
                this->_vtxInfo[vindex] = new typename Codim1Extractor<GV>::VertexInfo(vertex_index, vptr);
                // remember the vertex as a corner of the current face in temp_faces
                temp_faces.back().corners[i].idx = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only insert the index into the simplices array
                temp_faces.back().corners[i].idx = vimit->second->idx;
              }
            }

            // now increase the current face index
            simplex_index++;
            break;
          case 3 :
            assert(dim == 3);
            // we have a triangle here

            // add a new face to the temporary collection
            temp_faces.push_back(FaceInfo(simplex_index, eindex, *sit, FaceInfo::triangle));

            // try for each of the faces vertices whether it is already inserted or not
            for (int i = 0; i < simplex_corners; ++i)
            {
              // get the number of the vertex in the parent element
              int vertex_number = orientedSubface<dim>(gt, *sit, i);

              // get the vertex pointer and the index from the index set
              VertexPtr vptr(elit->template subEntity<dim>(vertex_number));
              IndexType vindex = this->indexSet().template index<dim>(*vptr);

              // remember the vertex' number in parent element's vertices
              temp_faces.back().corners[i].num = vertex_number;

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename Codim1Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
              if (vimit == this->_vtxInfo.end())
              {
                // insert into the map
                this->_vtxInfo[vindex] = new typename Codim1Extractor<GV>::VertexInfo(vertex_index, vptr);
                // remember the vertex as a corner of the current face in temp_faces
                temp_faces.back().corners[i].idx = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only insert the index into the simplices array
                temp_faces.back().corners[i].idx = vimit->second->idx;
              }
            }

            // now increase the current face index
            simplex_index++;
            break;
          case 4 :
            assert(dim == 3);
            // we have a quadrilateral here
            unsigned int vertex_indices[4];
            unsigned int vertex_numbers[4];

            // get the vertex pointers for the quadrilateral's corner vertices
            // and try for each of them whether it is already inserted or not
            for (int i = 0; i < cube_corners; ++i)
            {
              // get the number of the vertex in the parent element
              vertex_numbers[i] = orientedSubface<dim>(gt, *sit, i);

              // get the vertex pointer and the index from the index set
              VertexPtr vptr(elit->template subEntity<dim>(vertex_numbers[i]));
              IndexType vindex = this->indexSet().template index<dim>(*vptr);

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename Codim1Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
              if (vimit == this->_vtxInfo.end())
              {
                // insert into the map
                this->_vtxInfo[vindex] = new typename Codim1Extractor<GV>::VertexInfo(vertex_index, vptr);
                // remember this vertex' index
                vertex_indices[i] = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only remember the vertex' index
                vertex_indices[i] = vimit->second->idx;
              }
            }

            // now introduce the two triangles subdividing the quadrilateral
            // ATTENTION: the order of vertices given by "orientedSubface" corresponds to the order
            // of a Dune quadrilateral, i.e. the triangles are given by 0 1 2 and 3 2 1

            // add a new face to the temporary collection for the first tri
            temp_faces.push_back(FaceInfo(simplex_index++, eindex, *sit, FaceInfo::first));
            temp_faces.back().corners[0].idx = vertex_indices[0];
            temp_faces.back().corners[1].idx = vertex_indices[1];
            temp_faces.back().corners[2].idx = vertex_indices[2];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[0];
            temp_faces.back().corners[1].num = vertex_numbers[1];
            temp_faces.back().corners[2].num = vertex_numbers[2];

            // add a new face to the temporary collection for the second tri
            temp_faces.push_back(FaceInfo(simplex_index++, eindex, *sit, FaceInfo::second));
            temp_faces.back().corners[0].idx = vertex_indices[3];
            temp_faces.back().corners[1].idx = vertex_indices[2];
            temp_faces.back().corners[2].idx = vertex_indices[1];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[3];
            temp_faces.back().corners[1].num = vertex_numbers[2];
            temp_faces.back().corners[2].num = vertex_numbers[1];
            break;
          default :
            DUNE_THROW(Dune::NotImplemented, "the extractor does only work for triangle and quadrilateral faces (" << face_corners << " corners)");
            break;
          }
        }                         // end loop over found surface parts
      }
    }             // end loop over elements

    std::cout << "added " << simplex_index << " subfaces\n";

    // allocate the array for the face specific information...
    this->_faces.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());
  }


  // now first write the array with the coordinates...
  this->_coords.resize(this->_vtxInfo.size());
  typename Codim1Extractor<GV>::VertexInfoMap::const_iterator it1 = this->_vtxInfo.begin();
  for (; it1 != this->_vtxInfo.end(); ++it1)
  {
    // get a pointer to the associated info object
    typename Codim1Extractor<GV>::CoordinateInfo* current = &this->_coords[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }

}


template<typename GV>
inline void GeneralSurfaceExtractor<GV>::globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const
{
  // only interpolate barycentric in the given triangle (=> for flat quads this is exact, not so for non-flat quads!)
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  interpolateBarycentric<dimworld, ctype, Dune::FieldVector<ctype, dimworld> >(corners, bcoords, wcoords, dimworld);
}


template<typename GV>
inline void GeneralSurfaceExtractor<GV>::localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const
{
  if (this->_faces[index].category == 0)
  {
    // this is a triangle face
    Dune::array<Coords, simplex_corners> corners;
    unsigned int num_in_self = this->indexInInside(index);
    Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->type();
    for (int i = 0; i < simplex_corners; ++i)
      corners[i] = cornerLocalInRefElement<ctype, dimworld>(gt, num_in_self, i);
    interpolateBarycentric<dimworld, ctype, Dune::FieldVector<ctype, dimworld> >(corners, bcoords, ecoords, dimworld);
  }
  else
  {
    // this is a quadrilateral face
    Coords wcoords;
    this->localAndGlobalCoords(index, bcoords, ecoords, wcoords);
  }
}


template<typename GV>
inline void GeneralSurfaceExtractor<GV>::localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const
{
  this->globalCoords(index, bcoords, wcoords);
  // for rectangles avoid using world coordinates
  if (this->_faces[index].category == 0)
    // triangle face
    this->localCoords(index, bcoords, ecoords);
  else
    // arbitrarily distorted quadrilateral face (i.e. non-rectangular)
    ecoords = this->_elmtInfo.find(this->_faces[index].parent)->second->p->geometry().local(wcoords);
}


template<typename GV>
template<typename CoordContainer>
void GeneralSurfaceExtractor<GV>::globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords) const
{
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  for (int i = 0; i < wcoords.size(); ++i)
    interpolateBarycentric<dimworld, ctype, Dune::FieldVector<ctype, dimworld> >(corners, bcoords[i], wcoords[i], dimworld);

}


template<typename GV>
template<typename CoordContainer>
void GeneralSurfaceExtractor<GV>::localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords) const
{
  if (this->_faces[index].category == 0)
  {
    // this is a triangle face
    Dune::array<Coords, simplex_corners> corners;
    unsigned int num_in_self = this->indexInInside(index);
    Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->type();
    for (int i = 0; i < simplex_corners; ++i)
      corners[i] = cornerLocalInRefElement<ctype, dimworld>(gt, num_in_self, i);
    for (int i = 0; i < ecoords.size(); ++i)
      interpolateBarycentric<dimworld, ctype, Dune::FieldVector<ctype, dimworld> >(corners, bcoords[i], ecoords[i], dimworld);
  }
  else
  {
    // this is a quadrilateral face
    Dune::array<Dune::FieldVector<ctype,dim-1>, dimworld> subEntityCoords;
    for (int i=0; i<dimworld; i++)
      subEntityCoords[i] = barycentricToReference(bcoords[i]);

    CoordContainer wcoords;
    this->localAndGlobalCoords(index, subEntityCoords, ecoords, wcoords);
  }
}


template<typename GV>
void GeneralSurfaceExtractor<GV>::localAndGlobalCoords(unsigned int index,
                                                       const Dune::array<Dune::FieldVector<ctype,dim-1>, dimworld> &subEntityCoords,
                                                       Dune::array<Dune::FieldVector<ctype,dim>, dimworld> &ecoords,
                                                       Dune::array<Dune::FieldVector<ctype,dimworld>, dimworld> &wcoords) const
{
  Dune::array<Dune::FieldVector<ctype,dim>, dimworld> barycentricCoords;
  for (int i=0; i<dimworld; i++)
    barycentricCoords[i] = referenceToBarycentric(subEntityCoords[i]);

  this->globalCoords(index, barycentricCoords, wcoords);
  // for triangles
  if (this->_faces[index].category == 0)
    this->localCoords(index, barycentricCoords, ecoords);
  else
  {
    ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].parent)->second->p;
    for (int i = 0; i < ecoords.size(); ++i)
      ecoords[i] = eptr->geometry().local(wcoords[i]);
  }
}

#endif // GENERALSURFACEEXTRACTOR_HH_