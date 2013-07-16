#ifndef PHYSICS_ASSEMBLY_H
#define PHYSICS_ASSEMBLY_H


#include "agnosDefines.h"



#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/fe.h"

#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

namespace AGNOS
{


  /********************************************//**
   * \brief Base physics assembly class for deriving new physics classes.
   * Derived from libmesh Assembly class.
   ***********************************************/
  template<class T_S>
  class PhysicsAssembly : public libMesh::System::Assembly
  {
    public:
      PhysicsAssembly(
        libMesh::EquationSystems& es, 
        const std::string&        systemName
          );
      ~PhysicsAssembly( );

      void setParameterValues( const T_S& parameters );

      virtual void assemble( ) = 0;

      virtual void setSystemData( const GetPot& input ) = 0;

    protected:
      T_S                       m_paramValues;
      libMesh::EquationSystems& m_es; 
      const std::string&        m_systemName;

  };

  template<class T_S>
    PhysicsAssembly<T_S>::PhysicsAssembly( 
        libMesh::EquationSystems& es, 
        const std::string&        systemName
        )
    : m_es(es), m_systemName(systemName), m_paramValues( T_S(1) )
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    PhysicsAssembly<T_S>::~PhysicsAssembly( ){ }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    void PhysicsAssembly<T_S>::setParameterValues( const T_S& parameters) 
    { 
      m_paramValues = parameters ;
      return; 
    }
  

}


#endif // PHYSICS_ASSEMBLY_H
