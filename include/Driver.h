
#ifndef DRIVER_H
#define DRIVER_H


#include "Element.h"
#include <mpi.h>
#include <mpi.h>
#include <H5Cpp.h>



namespace AGNOS
{
  class H5IO;

  /********************************************//**
   * \brief Base class for driving the construction of SurrogateModel objects
   ***********************************************/
  class Driver
  { 

    // TODO comment variables and methods
    public:

      Driver( );
      Driver( const Communicator& comm, const Communicator& physicsComm, GetPot& input );

      virtual ~Driver( );

      void build( );
      void buildFromRestart( );
      void buildEvaluator( H5IO& h5io );

      void run( ) ;

      T_P evaluate( std::string, T_S& parameterValues ) ;

      void postProcess( ) ;
      void printSolutionData( ) ;
      void printSurrogateSettings(  ) ;
      void printParameterSettings(  ) ;
      void printDriverSettings(   ) ;
      void printSettings( ) ;
      void printSolution( unsigned int iteration=1 ) ;

      void computeErrors( 
          AGNOS::Element<T_S,T_P>&  elem,
          double& globalPhysics,
          double& globalTotal,
          double& globalSurrogate,
          double& maxElemError
          ) ;


    protected:
      GetPot& _input;
      std::shared_ptr< PhysicsModel<T_S,T_P> > _initPhysics( GetPot& input );
      std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > 
        _initSurrogate( GetPot& input, 
            std::vector< std::shared_ptr<AGNOS::Parameter> > parameters,
            std::shared_ptr<PhysicsModel<T_S,T_P> > physics) ;

      const Communicator& _comm;
      const Communicator& _physicsComm;

      // ---------------------
      // DRIVER VARIABLES
      /** maximum driver iterations */
      unsigned int  _maxIter;
      /** Determine which space to refine based on relative error estiamtes */
      bool          _adaptiveDriver;
      double        _refinePercentage ;
      /** force simultaneous refinement of both spaces */
      bool _simultRefine;
      /** restart HDF5 file */
      H5IO* _h5io;
      // ---------------------
      // ---------------------
      
      // ---------------------
      // ELEMENT VARIABLES
      std::list<AGNOS::Element<T_S,T_P> > _activeElems;
      std::queue<AGNOS::Element<T_S,T_P> >        _elemsToUpdate;
      // ---------------------
      

      // ---------------------
      // PARAMETERS VARIABLES
      unsigned int              _paramDim;
      std::vector< std::shared_ptr<AGNOS::Parameter> > _globalParameters;
      /** number of times to initially h refine the parameter domain before
       * doing any computation */
      unsigned int              _nInitialHRefinements ;
      // ---------------------
      
      // ---------------------
      // SURROGATE VARIABLES
      std::map<std::string, unsigned int> _surrogateNames;
      bool _refineSurrogate;
      bool _hRefine;
      bool _pRefine;
      std::vector<unsigned int> _pIncrement;
      bool _anisotropic;
      // ---------------------

      // ---------------------
      // PHYSICS VARIABLES
      bool _refinePhysics;
      bool _uniformRefine;
      // ---------------------
      

      // ---------------------
      // OUTPUT VARIABLES
      // TODO control on per surrogate basis
      std::string               _outputFilename; 
      std::shared_ptr<std::ostream>   _os ;
      std::vector<std::string>  _solutionsToPrint ;
      bool                      _computeMeans;
      bool                      _computeNorms;
      bool                      _outputIterations  ;
      bool                      _outputCoefficients  ;
      bool                      _outputErrorCoefficients  ;
      bool                      _outputWeights       ;
      bool                      _outputPoints        ;
      bool                      _outputIndexSet      ;

      bool                      _generateSamples      ;
      std::string               _sampleFile ;
      unsigned int              _nSamples ;

  };

  
}
#endif // DRIVER_H


