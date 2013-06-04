
#ifndef PHYSICS_VISCOUS_BURGERS_H
#define PHYSICS_VISCOUS_BURGERS_H

#include "agnosDefines.h"
#include "PhysicsModel.h"

namespace AGNOS
{
  /********************************************//**
   * \brief Basic 1D Burger's PhysicsModel class
   *
   * This example is given in the book
   * "Spectral Methods for Uncertainty Quantification" by Le Maitre and Knio
   ***********************************************/
  template<T_S, T_P>
  class PhysicsViscousBurgers : public PhysicsModel<T_S,T_P>
  {
    public:
      PhysicsViscousBurgers( GetPot& physicsInput ) ;
      
      ~PhysicsViscousBurgers( ) ;

    private:
      GetPot  m_input;
      double     m_min;
      double     m_max;
      double  m_viscosity;
  };


  template<T_S, T_P>
    PhysicsViscousBurgers::PhysicsViscousBurgers(
        GetPot& physicsInput
        )
      : m_input(physicsInput)
    {
      std::cout << "\nReading Physics model data\n";

      m_min = m_input("min",-10.);
      m_max = m_input("max",10.);
      m_viscosity = m_input("viscosity",1.);

      std::cout << "\n Initializing libmesh model\n"
      LibMeshInit init (argc, argv);

      return;
    }

  template<T_S, T_P>
    PhysicsViscousBurgers::~PhysicsViscousBurgers( )
    {
      return;
    }


  template<T_S, T_P>
    void PhysicsViscousBurgers::solvePrimal(
        T_S& paramaterValue)
    {
      // TODO   solve in libmesh the following problem
      //        find u \in U = {h1 s.t. u(-10)=0, u(10)=1} \approx ?
      //        a(u,v)=0 
      //        a(u,v) = int_{-10}^{10} [ u(1-u) - 2 \mu u_x ] v_x dx
      //        forall v \in V = {h1 s.t. v(-10)=v(10)=0}
      m_primalSolution = 1.0;
      return ;
    }


  template<T_S, T_P>
    void PhysicsViscousBurgers::solveAdjoint(
        T_S& paramaterValue,
        T_P& primalSolution
        )
    {
      // TODO   solve in libmesh the following problem
      //        find z \in V = {h1 s.t. z(-10)=0, z(10)=0}
      //        a'(u;v,z)=0 
      //        a'(u;v,z) = int_{-10}^{10} [ (1-2u)v - 2 \mu v_x ] z_x dx
      //        forall v \in U 
      m_adjointSolution = -1.0;
      return ;
    }



  template<T_S, T_P>
    double PhysicsViscousBurgers::evaluateQoi( 
        T_S& parameterValue,
        T_P& primalSolution    
        )
    {
      return 10.0;
    }

  template<T_S, T_P>
    double PhysicsViscousBurgers::estimateError( 
        T_S& parameterValue,  
        T_P& primalSolution,   
        T_P& adjointSolution  
        )
    {
      return 20.0;
    }

}

#endif // PHYSICS_VISCOUS_BURGERS_H
