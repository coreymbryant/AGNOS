
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
      PhysicsModel( 
          const Communicator& comm_in,
          const GetPot& input
          ) :
        _input(input),
        _communicator(comm_in),
        _compute_function(NULL)
      {
        // -----------------------------------------------------------
        // which solutions do we want to compute a surrogate for
        _solutionNames.clear( );

        for(unsigned int i=0; i<input.vector_variable_size("physics/solutions") ; i++)
          _solutionNames.insert( input("physics/solutions"," ",i) ) ;
        
        // default to only primal solution
        if(_solutionNames.size() == 0)
          _solutionNames.insert("primal");

        if(AGNOS_DEBUG)
        {
          std::set<std::string>::iterator it = _solutionNames.begin();
          std::cout << "Requested solutons: " ;
          for (; it!=_solutionNames.end(); ++it)
            std::cout <<  *it << " " ;
          std::cout << std::endl ;
        }


        if(AGNOS_DEBUG)
          std::cout << "solutionNames.size():" << _solutionNames.size() << std::endl;
        // -----------------------------------------------------------


      }

      /** Destructor */
      virtual ~PhysicsModel( ){};  /**< Default destructor */

      /** Function called by SurrogateModel to solve for requested solution
       * vectors at each evaluation point in parameter space  */
      virtual void compute( 
          const T_S& paramVector, 
          std::map<std::string, T_P >& solutionVectors 
          ) 
      { 
        if(AGNOS_DEBUG)
          std::cout << "DEBUG: calling provided compute function\n" ;

        if(_compute_function != NULL)
          this->_compute_function(paramVector, solutionVectors);

        if(AGNOS_DEBUG)
          std::cout << "DEBUG: ending call to provided compute function\n" ;

      }

      /** Refinement methods for the physics model */
      virtual void refine( ) {};

      /** return requested solution names */
      std::set<std::string> getSolutionNames( ) const
        { return _solutionNames; }


      /** attach a user defined compute function */
      void attach_compute_function(
          void fptr(const T_S& paramVector,std::map<std::string,T_P>& solutionVectors)
          )
      {
        libmesh_assert(fptr);
        _compute_function = fptr;
      }
      
    protected:
      /** communicator reference */
      const Communicator &_communicator;
      /** reference to input file just incase its needed after initialization */
      const GetPot&         _input;
      /** set of solutions to get (e.g. "primal","adjoint","qoi",etc. */
      std::set<std::string> _solutionNames;

      /** Function pointer to compute function */
      void(* _compute_function)(
          const T_S& paramVector, 
          std::map<std::string, T_P >& solutionVectors ) ;
      
  }; // PhysicsModel class

}

#endif //PHYSICS_MODEL_H
