
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

      /** Default constructor for abastract base class. Default to
       * solutionNames=["simple"] since we may not have primal, adjoint, qoi,
       * etc. */
      PhysicsModel( const Communicator& comm_in ) :
        _communicator(comm_in)
      {
        _solutionNames.insert("simple");
      }

      /** Destructor */
      virtual ~PhysicsModel( ){};  /**< Default destructor */

      /** Function called by SurrogateModel to solve for requested solution
       * vectors at each evaluation point in parameter space  */
      virtual void compute( 
          const T_S& paramVector, 
          std::map<std::string, T_P > solutionVectors 
          ) = 0;

      /** Refinement methods for the physics model */
      virtual void refine( ) = 0;
      std::set<std::string> getSolutionNames( ) const
        { return _solutionNames; }
      
    protected:
      /** communicator reference */
      const Communicator &_communicator;
      /** set of solutions to get (e.g. "primal","adjoint","qoi",etc. */
      std::set<std::string> _solutionNames;

      
  }; // PhysicsModel class

}

#endif //PHYSICS_MODEL_H
