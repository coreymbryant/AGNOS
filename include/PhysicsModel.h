

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
        _communicator(comm_in)
      {
        // -----------------------------------------------------------
        // which solutions do we want to compute a surrogate for
        _solutionNames.clear( );

        for(unsigned int i=0; i<input.vector_variable_size("solutions") ; i++)
          _solutionNames.insert( input("solutions"," ",i) ) ;
        
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
      virtual ~PhysicsModel( ){};  

      /** Function called by SurrogateModel to solve for requested solution
       * vectors at each evaluation point in parameter space. If solutionVectors
       * is non-empty then this routine should use provided soltuions in place
       * of any requested solutionNames. This allows for the construction of an
       * errorSurrogate where surrogate primal and adjoint solutions are
       * provided in solutionVectors.
       *
       * Must be redefined in derived classes*/
      virtual void compute( 
          const T_S& paramVector, std::map<std::string, T_P >& solutionVectors 
          ) = 0;

      /** Refinement methods for the physics model */
      virtual void refine( ) {};

      /** return requested solution names */
      std::set<std::string> getSolutionNames( ) const
        { return _solutionNames; }

      /** return reference to communicator */
      const Communicator& comm() const { return _communicator; }
      
      /** return reference to GetPot input object */
      const GetPot& input() const { return _input; }
      
    protected:
      /** communicator reference */
      const Communicator &_communicator;
      /** reference to input file just incase its needed after initialization */
      const GetPot&         _input;
      /** set of solutions to get (e.g. "primal","adjoint","qoi",etc. */
      std::set<std::string> _solutionNames;

      /** derived PhysicsModel classes need to handle settig parameter values
       * themselves */
      virtual void _setParameterValues( const T_S& parameterValues ) = 0;
      
  }; // PhysicsModel class

}

#endif //PHYSICS_MODEL_H
