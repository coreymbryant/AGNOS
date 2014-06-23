
#ifndef ELEMENT_H
#define ELEMENT_H

#include "PhysicsModel.h"
#include "Parameter.h"
#include "SurrogateModel.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Base element class
   *
   * Container for all necessary information/objects assocaited with an element
   * in stochastic space. Holds pointers to assocaited PhysiscModel,
   * SurrogateModel, and Parameters objects. 
   *
   * TODO: 
   *  load balancing
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class Element
  {

    public:
      /** Default constructor */
      Element( 
          std::vector< std::shared_ptr<AGNOS::Parameter> >&       parameters,
          std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > >&
          surrogates, 
          std::shared_ptr< PhysicsModel<T_S,T_P> >&               physics,
          double weight = 1.0
          ) ;

      /** Default destructor */
      virtual ~Element( ) ;

      /** determine if a parameter point lies in this element */
      bool covers( T_S& parameterValue )
      {
        bool result = true;
        for(unsigned int i=0;i<_parameters.size();i++)
        {
          if (_parameters[i]->type() == CONSTANT ) 
          {
            if ( parameterValue(i) != _parameters[i]->min() ) 
            {
              std::cout << "parameterValue(" << i << "): " 
                << parameterValue(i) << std::endl;
              std::cout << " param_min: " << _parameters[i]->min() << std::endl;
              std::cout << " param_max: " << _parameters[i]->max() << std::endl;
              result = false;
              break;
            }
          }
          else if (
              (parameterValue(i) < _parameters[i]->min())
              || (parameterValue(i) >= _parameters[i]->max()) 
             )
          {
            std::cout << "parameterValue(" << i << "): " 
              << parameterValue(i) << std::endl;
            std::cout << " param_min: " << _parameters[i]->min() << std::endl;
            std::cout << " param_max: " << _parameters[i]->max() << std::endl;
            result = false;
            break;
          }
        }

        return result;
      }

      /** split element into 2^dim children */
      std::vector< Element<T_S,T_P> > split( ) ;

      /** return parameters reference */
      std::vector<std::shared_ptr<AGNOS::Parameter> >   parameters( ) const
      { return _parameters; }

      /** return SurrogateModel pointer */
      std::vector< std::shared_ptr< SurrogateModelBase<T_S,T_P> > > surrogates(
          ) const { return _surrogates; }

      /** return PhysicsModel pointer */
      std::shared_ptr<PhysicsModel<T_S,T_P> >  physics( ) const
      { return _physics; }

      /** Set surrogate model to new model */
      void setSurrogates(
          std::vector< std::shared_ptr< SurrogateModelBase<T_S,T_P> > >&
          newSurrogates )
      { _surrogates = newSurrogates; }

      /** Set physics pointer to a different model */
      void setPhysics(
          std::shared_ptr< PhysicsModel<T_S,T_P> >& newPhysics )
      { _physics = newPhysics; }

      double _physicsError ;
      double _surrogateError ;
      double _totalError ;

      /** get weight of this element */
      double weight(){ return _weight;}

    protected:
      /** reference for parameters */
      std::vector< std::shared_ptr<AGNOS::Parameter> >  _parameters;

      /** pointer to SurrogateModels */
      std::vector< std::shared_ptr< SurrogateModelBase<T_S,T_P> > > _surrogates;

      /** pointer to PhysicsModel */
      std::shared_ptr< PhysicsModel<T_S,T_P> > _physics;

      /** weight of this element. Needed in multielement computations */
      double _weight;

  }; // Element class

  /********************************************//**
   * \brief utility function to compute means over a selection of elements
   ***********************************************/
  void computeMeans( 
      std::vector<std::string>&  solutions,
      std::list<AGNOS::Element<T_S,T_P> >& activeElems,
      std::map<std::string,T_P>& globalMeans 
      );
  void computeNorms( 
      std::vector<std::string>&  solutions,
      std::list<AGNOS::Element<T_S,T_P> >& activeElems,
      std::map<std::string,T_P>& globalNorms 
      );


}
#endif // ELEMENT_H
