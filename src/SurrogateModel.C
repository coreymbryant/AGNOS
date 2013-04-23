
#include "SurrogateModel.h"

namespace AGNOS
{

  SurrogateModel::SurrogateModel( )
  {
    return;
  }

  SurrogateModel::~SurrogateModel( )
  {
    return;
  }

  void SurrogateModel::setParameterDimension( unsigned int parameterDimension )
  {
    m_parameterDimension = parameterDimension ;
    return;
  }

  unsigned int SurrogateModel::getParameterDimension( ) const
  {
    return m_parameterDimension ;
  }

}
