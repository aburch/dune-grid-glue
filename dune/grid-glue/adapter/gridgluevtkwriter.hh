// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueVtkWriter.hh
 *  Version:     1.0
 *  Created on:  Mar 5, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Class thought to make graphical debugging of couplings easier.
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Write all remote intersections to a vtk file for debugging
 */

#ifndef GRIDGLUEVTKWRITER_HH_
#define GRIDGLUEVTKWRITER_HH_


#include <fstream>
#include <iomanip>
#include <vector>
#include <list>

#include <dune/common/geometrytype.hh>
#include <dune/grid/common/genericreferenceelements.hh>


const char TypeNames[][8] = { "double", "float  ", "int    " };

// double is default
template<typename T>
struct Nametraits
{
  enum { nameidx = 0 };
};


template<>
struct Nametraits<double>
{
  enum { nameidx = 0 };
};

template<>
struct Nametraits<float>
{
  enum { nameidx = 1 };
};

template<>
struct Nametraits<int>
{
  enum { nameidx = 2 };
};


/** \brief Write remote intersections to a vtk file for debugging purposes
 */
class GridGlueVtkWriter
{

  /** \brief Write either the grid0 or the grid1-side into streams
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
   */
  template <class Glue, int side>
  static void writeOneSide(const Glue& glue, const std::string& fngrid, const std::string& fnmerged)
  {
    dune_static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    std::ofstream fgrid, fmerged;

    fgrid.open(fngrid.c_str());
    fmerged.open(fnmerged.c_str());

    typedef typename Dune::SelectType<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::Type GridView;
    typedef typename GridView::Traits::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Dune::SelectType<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::Type RemoteIntersectionIterator;

    typedef typename GridView::ctype ctype;

    const int dim = GridView::dimension;
    const int dimw = Glue::dimworld;
    const int domdimw = GridView::dimensionworld;

    // coordinates have to be in R^3 in the VTK format
    std::string coordinatePadding;
    for (int i=domdimw; i<3; i++)
      coordinatePadding += " 0";

    const GridView& gridView = glue.template gridView<side>();

    int overlaps = glue.merger()->nSimplices();

    fgrid << "# vtk DataFile Version 2.0\nFilename: " << fngrid << "\nASCII" << std::endl;

    // WRITE POINTS
    // ----------------
    typedef typename Glue::Grid0Patch Extractor;
    std::vector<typename Extractor::Coords> coords;
    glue.template patch<side>().getCoords(coords);

    fgrid << "DATASET POLYDATA\nPOINTS " << coords.size() << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;

    for (size_t i=0; i<coords.size(); i++)
      fgrid << coords[i] << coordinatePadding << std::endl;

    // the merged grid (i.e. the set of remote intersections
    fmerged << "# vtk DataFile Version 2.0\nFilename: " << fnmerged << "\nASCII" << std::endl;
    fmerged << "DATASET POLYDATA\nPOINTS " << overlaps*(dim+1) << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;

    for (RemoteIntersectionIterator isIt = glue.template ibegin<side>();
         isIt != glue.template iend<side>();
         ++isIt)
    {
      for (int i = 0; i < isIt->geometry().corners(); ++i)
        fmerged << isIt->geometry().corner(i) << coordinatePadding << std::endl;
    }

    // WRITE POLYGONS
    // ----------------

    std::vector<typename Extractor::VertexVector> faces;
    std::vector<Dune::GeometryType> geometryTypes;
    glue.template patch<side>().getFaces(faces);
    glue.template patch<side>().getGeometryTypes(geometryTypes);

    unsigned int faceCornerCount = 0;
    for (size_t i=0; i<faces.size(); i++)
      faceCornerCount += faces[i].size();

    fgrid << "POLYGONS " << geometryTypes.size() << " " << geometryTypes.size() + faceCornerCount << std::endl;

    for (size_t i=0; i<faces.size(); i++) {

      fgrid << faces[i].size();

      // vtk expects the vertices to by cyclically ordered
      // therefore unfortunately we have to deal with several element types on a case-by-case basis
      if (geometryTypes[i].isSimplex()) {
        for (int j=0; j<dim; j++)
          fgrid << " " << faces[i][j];

      } else if (geometryTypes[i].isQuadrilateral()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2];

      } else {
        DUNE_THROW(Dune::NotImplemented, "Geometry type " << geometryTypes[i] << " not supported yet");
      }

      fgrid << std::endl;
    }

    fgrid << std::endl;

    int domainSimplexCorners = dim-Glue::Grid0Patch::codim+1;
    fmerged << "POLYGONS " << overlaps << " " << (domainSimplexCorners+1)*overlaps << std::endl;

    for (int i = 0; i < overlaps; ++i) {
      fmerged << domainSimplexCorners;
      for (int j=0; j<domainSimplexCorners; j++)
        fmerged << " " << domainSimplexCorners*i+j;
      fmerged << std::endl;
    }

#if 0
    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (gridSubEntityData.size()-1);

    fgrid << "CELL_DATA " << gridSubEntityData.size() << std::endl;
    fgrid << "SCALARS property_coding " << TypeNames[Nametraits<ctype>::nameidx] << " 1" << std::endl;
    fgrid << "LOOKUP_TABLE default" << std::endl;

    fmerged << "CELL_DATA " << overlaps << std::endl;
    fmerged << "SCALARS property_coding " << TypeNames[Nametraits<ctype>::nameidx] << " 1" << std::endl;
    fmerged << "LOOKUP_TABLE default" << std::endl;

    for (typename GridSubEntityData::const_iterator sEIt = gridSubEntityData.begin();
         sEIt != gridSubEntityData.end();
         ++sEIt, accum += delta)
    {
      // "encode" the parent with one color...
      fgrid << accum << std::endl;
      // ...and mark all of its merged grid parts with the same color
      for (int j = 0; j < sEIt->first.second; ++j)
        fmerged << accum << std::endl;
    }
#endif
    fgrid.close();
    fmerged.close();
  }

public:
  template<typename Glue>
  static void write(const Glue& glue, const std::string& filenameTrunk)
  {

    // Write extracted grid and remote intersection on the grid0-side
    writeOneSide<Glue,0>(glue,
                         filenameTrunk + "-domain.vtk",
                         filenameTrunk + "-intersections-domain.vtk");

    // Write extracted grid and remote intersection on the grid1-side
    writeOneSide<Glue,1>(glue,
                         filenameTrunk + "-target.vtk",
                         filenameTrunk + "-intersections-target.vtk");

  }

};



#endif // GRIDGLUEVTKWRITER_HH_
