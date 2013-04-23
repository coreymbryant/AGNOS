
#include "PhysicsModel.h"

namespace AGNOS
{
  /********************************************//**
   * \brief Basic 1D Burger's equation
   *
   * This example is given in the Le Maitre book
   ***********************************************/
  class ViscousBurgersPhysics : public PhysicsModel
  {
    public:
      ViscousBurgersPhysics( ) ;
      ViscousBurgersPhysics(int xMin, int xMin) ;
      
      ~ViscousBurgersPhysics( ) ;

    private:
      int xMin;
      int xMax;
      double m_viscosity;


  }
}
