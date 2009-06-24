// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUECOMMUNICATE_HH
#define DUNE_GRIDGLUECOMMUNICATE_HH

/**@file
   @author Christian Engwer
   @brief Describes the parallel communication interface class for Dune::GridGlue
 */

#include <dune/grid/common/datahandleif.hh>
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /**
     \brief describes the features of a data handle for
     communication in parallel runs using the GridGlue::communicate methods.
     Here the Barton-Nackman trick is used to interprete data handle objects
     as its interface. Therefore usable data handle classes need to be
     derived from this class.

     \tparam DataHandleImp implementation of the users data handle
     \tparam DataTypeImp type of data that are going to be communicated which is exported as \c DataType (for example double)
     \ingroup GICollectiveCommunication
   */
  template <class DataHandleImp, class DataTypeImp>
  class GridGlueCommDataHandleIF : public CommDataHandleIF<DataHandleImp, DataTypeImp>
  {
  protected:
    // one should not create an explicit instance of this inteface object
    GridGlueCommDataHandleIF() {}

  public:

    /** @brief pack data from user to message buffer
        @param buff message buffer provided by the grid
        @param e entity for which date should be packed to buffer
        @param i RemoteIntersection for which date should be packed to buffer
     */
    template<class MessageBufferImp, class EntityType, class RISType>
    void gather (MessageBufferImp& buff, const EntityType& e, const RISType & i) const
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().gather(buffIF,e)));
    }

    /*! unpack data from message buffer to user
       n is the number of objects sent by the sender
       @param buff message buffer provided by the grid
       @param e entity for which date should be unpacked from buffer
       @param n number of data written to buffer for this entity before
     */
    template<class MessageBufferImp, class EntityType, class RISType>
    void scatter (MessageBufferImp& buff, const EntityType& e, const RISType & i, size_t n)
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().scatter(buffIF,e,n)));
    }

  private:
    //!  Barton-Nackman trick
    DataHandleImp& asImp () {
      return static_cast<DataHandleImp &> (*this);
    }
    //!  Barton-Nackman trick
    const DataHandleImp& asImp () const
    {
      return static_cast<const DataHandleImp &>(*this);
    }
  }; // end class CommDataHandleIF

} // end namespace Dune
#endif