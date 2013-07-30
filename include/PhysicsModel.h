
// Do we need both versions of the defined functions, or just get by with the
// one given the solution and pass it in if needed
//
// Solve functions intentially do not return the solution. There may be cases
// where we want to solve the problem internally but not create a new instance
// of the solution. Create get..Solution functions to check if solution exists
// and if not solve problem first. These should be what is used by surrogate
// model.



#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H

#include "agnosDefines.h"
// libmesh includes
#include "libmesh/parallel.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Base physics model class
   *
   * Abstract framework for PhysicsModel classes. 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsModel
  {
    
    public: 

      PhysicsModel( const Parallel::Communicator& comm_in ) :
        _communicator(comm_in)
      {
        _solutionNames.insert("simple");
      }

      virtual ~PhysicsModel( ){};  /**< Default destructor */

      virtual void compute( 
          const T_S& paramVector, 
          std::map<std::string, T_P > solutionVectors 
          ) = 0;

      virtual void refine( ) = 0;
      std::set<std::string> getSolutionNames( ) const
        { return _solutionNames; }
      
    protected:
      const Parallel::Communicator &_communicator;
      std::set<std::string> _solutionNames;

      
  }; // PhysicsModel class

}

#endif //PHYSICS_MODEL_H
