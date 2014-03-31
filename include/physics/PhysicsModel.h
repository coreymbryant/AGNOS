

#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H

#include "agnosDefines.h"

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
        // which solutions are available in this phyisc class?
        // default to only primal solution
        if(_availableSolutions.size() == 0)
          _availableSolutions.insert("primal");

        if(AGNOS_DEBUG)
        {
          std::set<std::string>::iterator it = _availableSolutions.begin();
          std::cout << "Requested solutons: " ;
          for (; it!=_availableSolutions.end(); ++it)
            std::cout <<  *it << " " ;
          std::cout << std::endl ;
        }


        if(AGNOS_DEBUG)
          std::cout << "solutionNames.size():" << _availableSolutions.size() << std::endl;
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
          std::set<std::string>& computeSolutions,
          const T_S& paramVector, std::map<std::string, T_P >& solutionVectors 
          ) = 0;

      /** Refinement methods for the physics model */
      virtual void refine( ) {};

      /** return requested solution names */
      std::set<std::string> getAvailableSolutions( ) const
        { return _availableSolutions; }

      /** return reference to communicator */
      const Communicator& comm() const { return _communicator; }
      
      /** return reference to GetPot input object */
      const GetPot& input() const { return _input; }
      
    protected:
      /** communicator reference */
      const Communicator &_communicator;
      /** reference to input file just incase its needed after initialization */
      const GetPot&         _input;
      /** set of solutions available for this physics 
       * (e.g.  "primal","adjoint","qoi",etc. */
      std::set<std::string> _availableSolutions;

      /** derived PhysicsModel classes need to handle settig parameter values
       * themselves */
      virtual void _setParameterValues( const T_S& parameterValues ) = 0;
      
  }; // PhysicsModel class

}

#endif //PHYSICS_MODEL_H
