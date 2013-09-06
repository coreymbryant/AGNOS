
#ifndef ELEMENT_H
#define ELEMENT_H

/* #include "agnosDefines.h" */
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
          const Communicator&                                     comm,
          std::vector< std::shared_ptr<AGNOS::Parameter> >&       parameters,
          std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > >& surrogates, 
          std::shared_ptr< PhysicsModel<T_S,T_P> >                physics
          ) ;

      /** Default destructor */
      virtual ~Element( ) ;

      /** split element into 2^dim children */
      std::vector< Element<T_S,T_P> > split( ) ;

      // TODO do we need both or should we set a flag like in driver to control
      // what type of refinement is performed?
      // Do we even need these at all?
      /** refine element's surrogate model */
      /* void refinePhysics( ) ; */

      /** refine element's physics model */
      /* void refineSurrogate( ) ; */

      /** return comm reference */
      const Communicator& comm( ) const
      { return _comm; }

      /** return parameters reference */
      std::vector<std::shared_ptr<AGNOS::Parameter> >   parameters( ) const
      { return _parameters; }

      /** return SurrogateModel pointer */
      std::vector< std::shared_ptr< SurrogateModel<T_S,T_P> > > surrogates( ) const
      { return _surrogates; }

      /** return PhysicsModel pointer */
      std::shared_ptr<PhysicsModel<T_S,T_P> >  physics( ) const
      { return _physics; }


    private:
      /** reference to communicator */
      const Communicator& _comm;

      /** reference for parameters */
      std::vector< std::shared_ptr<AGNOS::Parameter> >  _parameters;

      /** pointer to SurrogateModels */
      std::vector< std::shared_ptr< SurrogateModel<T_S,T_P> > > _surrogates;

      /** pointer to PhysicsModel */
      std::shared_ptr< PhysicsModel<T_S,T_P> > _physics;

  }; // Element class


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
    Element<T_S,T_P>::Element(
        const Communicator&                                     comm,
        std::vector< std::shared_ptr<AGNOS::Parameter> >&       parameters,
        std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > >& surrogates, 
        std::shared_ptr< PhysicsModel<T_S,T_P> >                  physics
        ) 
    :
      _comm(comm),
      _parameters(parameters),
      _surrogates(surrogates),
      _physics(physics)
    {
    }
  
  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S, class T_P>
    Element<T_S,T_P>::~Element( )
    {
    }

  /********************************************//**
   * \brief Split element into new elements by bisection in each parameter
   * direction.
   *
   * Note: SurrogateModel parameters will not be updated, its up to the user to
   * decide what will be done after splitting the element, i.e. updating the
   * SurrogateModel or reuse old model and just restrict evaluation to new
   * parameter domain. 
   ***********************************************/
  template<class T_S, class T_P>
    std::vector< Element<T_S,T_P> > Element<T_S,T_P>::split( )
    {
      unsigned int dim = _parameters.size();
      unsigned int nChildren = std::pow(2,dim) ;

      // TODO 
      //  Chnage to vector of smart ptrs instead?
      // initialize vector for new elements
      std::vector< Element<T_S,T_P> > newElements;
      newElements.reserve( nChildren ) ;

      // TODO 
      //  split parameters
      std::vector<std::shared_ptr<AGNOS::Parameter> > newParameters =
        _parameters;
      /* newParameters.reserve( dim ); */
      /* for (unsigned int i=0; i<dim; i++ ) */

      //  initalize new Elements with same pointers
      for (unsigned int i=0; i<nChildren; i++)
        newElements.push_back( 
            Element<T_S,T_P>(
              _comm,
              newParameters,
              _surrogates,
              _physics
              )
            );

      return newElements;
    }
  
  

}
#endif // ELEMENT_H
