
#include "Element.h"

namespace AGNOS
{


  /********************************************//**
   * \brief Default constructor 
   *
   * Note that this method assumes provided parameters represent global
   * parameter space, i.e. the element weight is the 1. 
   ***********************************************/
  template<class T_S, class T_P>
    Element<T_S,T_P>::Element(
        std::vector< std::shared_ptr<AGNOS::Parameter> >&       parameters,
        std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > >& surrogates, 
        std::shared_ptr< PhysicsModel<T_S,T_P> >&                 physics,
        double weight
        ) 
    :
      _parameters(parameters),
      _surrogates(surrogates),
      _physics(physics),
      _weight(weight)
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
      double weight = _weight;

      // TODO 
      //  Chnage to vector of smart ptrs instead?
      // initialize vector for new elements
      std::vector< Element<T_S,T_P> > newElements;
      newElements.reserve( nChildren ) ;

      /* std::vector<std::shared_ptr<AGNOS::Parameter> > newParameters = */
      /*   _parameters; */

      //  split parameters
      //  Create new 1d parameters that will be needed in tensor product
      std::vector< std::vector< std::shared_ptr<AGNOS::Parameter> > > allParams;
      allParams.reserve( nChildren ) ;
      // first direction (to initialize everything)
      {
        double min = _parameters[0]->min() ;
        double max = _parameters[0]->max() ;
        double midpoint = (min + max)/ 2.0;

        // adjust child weight: divide by parent size
        /* weight /= _parameters[0]->measure(min,max) ; */
        // multiply by child size
        weight *= _parameters[0]->measure(min,midpoint);

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
        
        // adjust child weight: divide by parent size
        /* weight /= _parameters[i]->measure(min,max) ; */
        // multiply by child size
        weight *= _parameters[i]->measure(min,midpoint);

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
              _physics,
              weight
              )
            );
        
      }


      return newElements;
    }

  /********************************************//**
   * \brief utility function to compute means over a selection of elements
   ***********************************************/
  void computeMeans( 
      std::vector<std::string>&  solutionNames,
      std::list<AGNOS::Element<T_S,T_P> >& activeElems,
      std::map<std::string,T_P>& globalMeans 
      )
  {

    globalMeans.clear();

    // loop through elements
    std::list<AGNOS::Element<T_S,T_P> >::iterator elit =
      activeElems.begin();
    for (; elit!=activeElems.end(); elit++)
    {
      // retrieve element mean coefficients
      std::map<std::string,T_P> elementMeans
        = elit->surrogates()[0]->mean( );

      // add contributions from this element to global means
      for(unsigned int i=0 ; i<solutionNames.size();i++)
      {
        agnos_assert( (elementMeans.count(solutionNames[i]) > 0) ) ;

        // scale by elem weight
        elementMeans[solutionNames[i]].scale( elit->weight() );
        
        if(globalMeans.count(solutionNames[i])==0)
          globalMeans.insert( 
              std::pair<std::string,T_P>(
                solutionNames[i],elementMeans[solutionNames[i]] )
              );
        else
          globalMeans[solutionNames[i]] += elementMeans[solutionNames[i]] ;


      }
        
    }

    return ;
  }

  template class
    Element<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}
