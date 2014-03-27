

#ifndef PHYSICS_CATENARY_H
#define PHYSICS_CATENARY_H

#include "PhysicsUser.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsUser class - catenary chain
   *
   * A simple, single dof, system useful for testing purposes
   *
   * 1 dof corresponds to middle node in two (linear) element approximation on
   * domain [0,1]
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenary : public PhysicsUser<T_S,T_P>
  {

  public:
    /** Default constructor Default */
    PhysicsCatenary( 
        const Communicator& comm_in,
        const GetPot& input
        );

    /** Destructor */
    ~PhysicsCatenary( );

    /** Function called by SurrogateModel to solve for requested solution
     * vectors at each evaluation point in parameter space  */
    virtual void compute( 
        std::set<std::string>& computeSolutions,
        const T_S& paramVector, std::map<std::string, T_P >& solutionVectors 
        ) ;

    /** Refinement methods for the physics model */
    virtual void refine( );
  
  protected:
    /** communicator reference */
    const Communicator &_communicator;
    /** reference to input file just incase its needed after initialization */
    const GetPot&         _input;

    /** forcing magnitude */
    double m_forcing;

  };




/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenary<T_S,T_P>::PhysicsCatenary(
        const Communicator& comm_in,
        const GetPot& input
      ) :
    PhysicsUser<T_S,T_P>(comm_in,input),
    _communicator(comm_in),_input(input)
  {
    m_forcing = input("forcing",-10.0) ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenary<T_S,T_P>::~PhysicsCatenary( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void PhysicsCatenary<T_S,T_P>::compute(
        std::set<std::string>& computeSolutions,
        const T_S& paramVector,
        std::map<std::string, T_P >& solutionVectors 
        ) 
    {
      if(AGNOS_DEBUG)
        std::cout << "entering compute routine" << std::endl;
      solutionVectors.clear();

      if (!solutionVectors.count("primal"))
      {
        T_P primalSolution(1);
        primalSolution(0) =  m_forcing / (8. *  paramVector(0) ) ;
        solutionVectors.insert( 
            std::pair<std::string,T_P>( "primal", primalSolution )
              );
        if(AGNOS_DEBUG)
          std::cout << "compute -- primal solution size:" << primalSolution.size()
            << std::endl;
      }

      if( computeSolutions.count("adjoint" ) 
          || computeSolutions.count("errorEstimator") 
          || computeSolutions.count("errorIndicators") 
          )
      {
        if ( solutionVectors.count("adjoint")  )
        {
        }
        else
        {
          T_P adjointSolution(1);
          adjointSolution(0) =  1.0 / (4. * paramVector(0))  ;

          if ( computeSolutions.count("adjoint") )
          {
            solutionVectors.insert( 
                std::pair<std::string,T_P>( "adjoint", adjointSolution )
                  );
          }
        }
        
        // ----------------------------------
        // ERROR ESTIMATE AND ERROR INDICATORS
        // Compute global error estimate if requested
        if ( computeSolutions.count("errorEstimate") ) 
        {
          // in this case the FE solution interpolates at x=1/2 so the QoI is
          // evaluated exactly
          T_P error(1);
          error.zero();
          solutionVectors.insert( 
                    std::pair<std::string,T_P>( "errorEstimate", error )
                      );
        }
        
        // compute indicators if they were requested
        if (  computeSolutions.count("errorIndicators") )
        {
        }

      } // end if adjoint|errorEstimate|errorIndicators
      
      // QUANTITY OF INTEREST
      if ( computeSolutions.count("qoi") )
      {
        T_P qoiValue(1);
        qoiValue(0) = solutionVectors["primal"](0) ;
        solutionVectors.insert( 
                  std::pair<std::string,T_P>( "qoi", qoiValue )
                    );
      }


      if(AGNOS_DEBUG)
      {
        std::cout << "solutionVectors['primal'].size():" 
          << solutionVectors["primal"].size() << std::endl;

        std::cout << "leaving compute routine" << std::endl;
      }
    }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P> 
    void PhysicsCatenary<T_S,T_P>::refine( ) 
    {
      return;
    }


}

#endif // PHYSICS_CATENARY_H
