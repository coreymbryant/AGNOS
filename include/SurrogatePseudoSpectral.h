
#ifndef SURROGATE_PSEUDO_SPECTRAL_H
#define SURROGATE_PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"
#include "QuadratureTensorProduct.h"
#include "petscmat.h"


namespace AGNOS
{

/********************************************//**
 * \brief SurrogateModel based on Pseudo-spectral projection. 
 *
 * This class provides the framework for constructing surrogate models using
 * non-intrusive spectral projection methods. 
 *
 * As of now only the derived class PseudoSpectralTensorProduct is operational
 * but it could be extended to other methods as well. Both isotropic and
 * non-isotropic polynomial orders are supported. 
 ***********************************************/
  template<class T_S, class T_P>
    class SurrogatePseudoSpectral : public SurrogateModel<T_S,T_P>
  {

    public:

      /** Constructor */
      SurrogatePseudoSpectral( 
          const Communicator&               comm,
          std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
          const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
          const std::vector<unsigned int>&  order,
          std::set<std::string> computeSolutions = std::set<std::string>()
          );


      /** Secondary Constructor. 
       *  ***DO NOT MISTAKE FOR A COPY CONSTRUCTOR ***
       *  Intended use is for constructing a secondary surrogate model using the
       *  primary model as an evaluating object in the build routine. 
       *  If additional inputs are defined it will
       * construct a new surrogate increasing the order and using
       * primarySurrogate to perform evaluations in the constructions */
      SurrogatePseudoSpectral( 
          std::shared_ptr<SurrogateModel<T_S,T_P> > primarySurrogate, 
          std::vector<unsigned int> increaseOrder = std::vector<unsigned int>(),
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::set<std::string>(),
          std::set<std::string> computeSolutions = std::set<std::string>()
          );

      virtual ~SurrogatePseudoSpectral( );

      void build( ) ;
      std::map< std::string, T_P > computeContribution( 
          std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
          unsigned int index
          );

      using SurrogateModel<T_S,T_P>::evaluate; 
      std::map<std::string, T_P> evaluate( 
          std::set< std::string >  solutionNames,
          T_S&                        parameterValues,
          bool saveLocal = true /**< save solution locally after evaluation*/
          ) const ;

      using SurrogateModel<T_S,T_P>::l2Norm;
      std::map< std::string, T_P> l2Norm(
        std::set< std::string > solutionNames
        );
      double l2NormDifference(
          SurrogateModel<T_S,T_P>& comparisonModel,
          std::string solutionName );

      // Manipulators
      unsigned int              getNIntegrationPoints( ) const;
      std::vector<T_S>          getIntegrationPoints( ) const;
      std::vector<double>       getIntegrationWeights( ) const;
      /** Return index set for this surrogate */
      const std::vector< std::vector< unsigned int> > 
                                indexSet( ) const;
      std::vector<double>       evaluateBasis( 
          const std::vector< std::vector<unsigned int> >& indexSet,
          T_S& parameterValues ) const;

      void                      printIntegrationWeights( std::ostream& out ) const;
      void                      printIntegrationPoints( std::ostream& out ) const;
      void                      printIndexSet( std::ostream& out ) const;
      void                      prettyPrintIntegrationWeights( ) const;
      void                      prettyPrintIntegrationPoints( ) const;
      void                      prettyPrintIndexSet( ) const;


    protected:

      unsigned int              _nIntegrationPoints;
      std::vector<T_S>          _integrationPoints ;
      std::vector<double>       _integrationWeights ;
      std::vector<unsigned int> _integrationIndices;

      std::vector< std::vector<unsigned int> > _indexSet;



  };



  
}
#endif // SURROGATE_PSEUDO_SPECTRAL_H


