#include "PhysicsModel.h"

namespace AGNOS
{

  PhysicsModel::PhysicsModel( )
  {
    return;
  }

  PhysicsModel::~PhysicsModel( )
  {
    return;
  }


  const PhysicsDataType& PhysicsModel::getPrimalSolution( ) const
  {
    return m_primalSolution;
  }

  const PhysicsDataType& PhysicsModel::getAdjointSolution( ) const
  {
    return m_adjointSolution;
  }

  const PhysicsDataType& PhysicsModel::getErrorIndicators( ) const
  {
    return m_errorIndicators;
  }

}
