
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
          std::vector< std::shared_ptr<AGNOS::Parameter> >&       parameters,
          std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > >& surrogates, 
          std::shared_ptr< PhysicsModel<T_S,T_P> >&               physics
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

      /** return parameters reference */
      std::vector<std::shared_ptr<AGNOS::Parameter> >   parameters( ) const
      { return _parameters; }

      /** return SurrogateModel pointer */
      std::vector< std::shared_ptr< SurrogateModel<T_S,T_P> > > surrogates( ) const
      { return _surrogates; }

      /** return PhysicsModel pointer */
      std::shared_ptr<PhysicsModel<T_S,T_P> >  physics( ) const
      { return _physics; }

      /** Set surrogate model to new model */
      void setSurrogates(
          std::vector< std::shared_ptr< SurrogateModel<T_S,T_P> > >&
          newSurrogates )
      { _surrogates = newSurrogates; }

      /** Set physics pointer to a different model */
      void setPhysics(
          std::shared_ptr< PhysicsModel<T_S,T_P> >& newPhysics )
      { _physics = newPhysics; }


    protected:
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
        std::vector< std::shared_ptr<AGNOS::Parameter> >&       parameters,
        std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > >& surrogates, 
        std::shared_ptr< PhysicsModel<T_S,T_P> >&                 physics
        ) 
    :
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

      /* std::vector<std::shared_ptr<AGNOS::Parameter> > newParameters = */
      /*   _parameters; */

      // TODO 
      //  split parameters
      //  Create new 1d parameters that will be needed in tensor product
      std::vector< std::vector< std::shared_ptr<AGNOS::Parameter> > > allParams;
      allParams.reserve( nChildren ) ;
      // first direction (to initialize everything)
      {
        double min = _parameters[0]->min() ;
        double max = _parameters[0]->max() ;
        double midpoint = (min + max)/ 2.0;


        allParams.push_back( 
            std::vector< std::shared_ptr<AGNOS::Parameter> >(
              1,
              std::shared_ptr<AGNOS::Parameter>(
                new AGNOS::Parameter( _parameters[0]->type(), min, midpoint) )
              )
            );
        allParams.push_back( 
            std::vector< std::shared_ptr<AGNOS::Parameter> >(
              1,
              std::shared_ptr<AGNOS::Parameter>(
                new AGNOS::Parameter( _parameters[0]->type(), midpoint, max) )
              )
            );
      }
      // add aditional directions by building on first direction
      for (unsigned int i=1; i<dim; i++ )
      {
        double min = _parameters[i]->min() ;
        double max = _parameters[i]->max() ;
        double midpoint = (min + max)/ 2.0;

        const unsigned int oldSize = allParams.size();
        for (unsigned int j=0; j<oldSize; j++)
        {
          // push back a copy of this row before we alter it 
          allParams.push_back( allParams[j] );
          
          // push back this directions lesser parameter to this row
          allParams[j].push_back( 
                std::shared_ptr<AGNOS::Parameter>(
                  new AGNOS::Parameter( _parameters[i]->type(), min, midpoint) )
              );

          // then add this directions greater param range to last row
          // (the one we added at the beginning of the loop)
          allParams.back().push_back( 
                std::shared_ptr<AGNOS::Parameter>(
                  new AGNOS::Parameter( _parameters[i]->type(), midpoint, max) )
              );
        }
      }

      // create new elements from new parameters
      for (unsigned int c=0; c<nChildren; c++)
      {
        newElements.push_back( 
            Element<T_S,T_P>(
              allParams[c],
              _surrogates,
              _physics
              )
            );
        
      }


      return newElements;
    }
  
  

}
#endif // ELEMENT_H
