
#ifndef ADAPTIVE_DRIVER_H
#define ADAPTIVE_DRIVER_H

#include "Driver.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Driver routine for basic adaptive surrogate construction
   *
   * This class handles the computation of the total, physical, and surrogate
   * error contributions based on the physics and surrogate models provided. It
   * determines which approximation to refine and calls the appropriate
   * refinement function of that model class. 
   *
   * Particular types of refinement can be controlled based on SurrogateModel
   * and PhysicsModel class's refine() function
   * 
   * TODO derived classes could control type of adaptivity or leave that up to
   * how surrogate/physics models are defined, i.e. 
   * uniform x uniform h
   * unifrom x adaptive h
   * uniform x uniform p
   * h       x uniform h
   * p       x uniform
   * .
   * .
   * .
   ***********************************************/
  class DriverAdaptive : public Driver
  {

    public:

      DriverAdaptive( );           /**< Default constructor */
      DriverAdaptive( GetPot& input );

      ~DriverAdaptive( );  /**< Default destructor */

      void          setMaxIter(unsigned int maxIter);
      unsigned int  getMaxIter( ) const;

      void    setErrorTolerance(double errorTolerance);
      double  getErrorTolerance( ) const;

      void  setPhysicsModel(PhysicsModel* physicsModel) ;
      const PhysicsModel* getPhysicsModel( ) const;

      void  setSurrogateModel(SurrogateModel* surrogateModel) ;
      const SurrogateModel* getSurrogateModel( ) const;


    protected:

      unsigned int  m_maxIter;
      double        m_errorTolerance;
      bool          m_refinePhysical;
      bool          m_refineSurrogate;

      double m_totalError;
      double m_physicalError;
      double m_surrogateError;
  
  }




}

#endif //ADAPTIVE_DRIVER_H
