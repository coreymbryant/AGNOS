/* #include "agnosDefines.h" */

#include "Driver.h"
#include "ioHandler.h"
#include "Element.h"
#include "SurrogateModel.h"
#include "PseudoSpectralTensorProduct.h"
#include "PhysicsModel.h"
#include "PhysicsViscousBurgers.h"
#include "PhysicsCatenary.h"
#include "PhysicsCatenaryLibmesh.h"
#include "PhysicsDiffusion.h"
#ifdef AGNOS_ENABLE_CHANNELFLOW
#include "PhysicsChannelFlow.h"
#endif // AGNOS_ENABLE_CHANNEL_FLOW
#ifdef AGNOS_ENABLE_GRINS
#include "PhysicsGrins.h"
#endif // AGNOS_ENABLE_GRINS
#include "PhysicsLibmesh.h"
#include <mpi.h>
#include <gsl/gsl_rng.h>

namespace AGNOS
{



/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver( 
      const Communicator& comm,
      const Communicator& physicsComm,
      GetPot& input 
      ) :
    _comm(comm), _physicsComm(physicsComm), _input(input) 
  {
    if(AGNOS_DEBUG)
      std::cout << "rank: " << _comm.rank() <<  std::endl;

    
    // DRIVER SETTINGS
    _maxIter = _input("driver/maxIter",1);
    /** assumes one secondary surrogate used for total error estimation */
    _adaptiveDriver = _input("driver/adaptive",false);
    std::cout << "adapt? " << _adaptiveDriver << std::endl;
    _refinePercentage = _input("driver/refinePercentage",0.20);
    
    // ADAPTIVE SETTINGS
    _simultRefine = _input("driver/simultRefine",false);
    _forceSurrogateRefine = _input("driver/forceSurrogateRefine",false);
    _forcePhysicsRefine = _input("driver/forcePhysicsRefine",false);

    // restart capabilities
    bool restart = _input("driver/restart",false);
    bool evaluator = _input("driver/evaluator",false);
    if (restart)
    {
      std::string restartFile = _input("restart/fileName","agnos.h5");
      _h5io = new H5IO( restartFile, H5F_ACC_RDONLY );
      buildFromRestart( );
    }
    else if (evaluator)
    {
      std::string evaluatorFile 
        = _input("driver/evaluatorFile","agnos.h5");
      _h5io = new H5IO( evaluatorFile, H5F_ACC_RDONLY );
      buildEvaluator( *_h5io );
    }
    else 
    {
      _h5io = new H5IO( "agnos.h5", H5F_ACC_TRUNC );
      build();
    }
    
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  Driver::~Driver( )
  {
    delete _h5io ;
  }


/********************************************//**
 * \brief build simulation driver based on input file
 * 
 ***********************************************/
 void  Driver::build( )
  {
    
    // PARAMETER SETTINGS
    _paramDim = _input("parameters/dimension", 1);
    _nInitialHRefinements = _input("parameters/nInitialHRefinements", 0);

    // read in parameer dimension and bounds
    // TODO if less than dim but greater than one -> error
    //      if one and dim > 1 set all as same
    _globalParameters.reserve(_paramDim);


    // construct parameter structures
    std::string paramType = "CONSTANT";
    double paramMin = 1.0;
    double paramMax = 1.0;
    // read in parameter type and bounds
    // if all directions are not provided we assume they match the last
    // provided directions sections
    for (unsigned int i=0; i < _paramDim; i++)
    {
      paramType = _input("parameters/types",paramType,i);
      paramMin  = _input("parameters/mins",paramMin,i); 
      paramMax  = _input("parameters/maxs",paramMax,i);
      _globalParameters.push_back( std::shared_ptr<AGNOS::Parameter>(
            new AGNOS::Parameter( paramType, paramMin, paramMax )) 
          );
    }

    
    // PHYSICS SETTINGS
    std::shared_ptr< PhysicsModel<T_S,T_P> > physics = _initPhysics( _input );


    // SURROGATE MODEL(s) SETTINGS
    _input.set_prefix("surrogateModels/") ;
    _refineSurrogate = _input("refine",false);
    _hRefine = _input("hRefine",false);
    _pRefine = _input("pRefine",true);
    // pIncrement vector (default to 1 in all directions)
    int pIncrementDim = _input.vector_variable_size("pIncrement");
    if (pIncrementDim == 1)
      for (unsigned int i=0; i < _paramDim; i++)
        _pIncrement.push_back( _input("pIncrement", 1) ) ;
    else
      for (unsigned int i=0; i < _paramDim; i++)
        _pIncrement.push_back( _input("pIncrement", 1, i) ) ;
    _anisotropic = _input("anisotropic",false);

    /* initialize surrogate model container */
    std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates =
      _initSurrogate( _input, _globalParameters, physics ); 


    // Construct a single initial element and add it to the update queue
    AGNOS::Element<T_S,T_P> baseElement(
          _globalParameters,
          surrogates,
          physics
          ) ;
    _elemsToUpdate.push(baseElement);


    //perform initial h refinements
    for (unsigned int i=0; i<_nInitialHRefinements; i++)
    {
      unsigned int nElems = _elemsToUpdate.size() ;
      for (unsigned int j=0; j < nElems; j++)
      {
        std::vector< Element<T_S,T_P> > children = _elemsToUpdate.front().split() ;
        _elemsToUpdate.pop();

        for (unsigned int c=0; c<children.size(); c++)
        {
          // cosntruct a new surrogate for this element
          std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > childSurrogates =
              _initSurrogate( _input, children[c].parameters(), children[c].physics() ) ;
          children[c].setSurrogates( childSurrogates ); 

          // save element to update queue
          _elemsToUpdate.push(children[c]) ;
        } // children loop
      } // nElems loop
    } // nInitialHRefinements loop

    

    // OUTPUT DATA SETTINGS
    _input.set_prefix("output/") ;
    _outputFilename      = _input("filename","cout");
    if ( _outputFilename == "cout" )
      _os.reset( &(std::cout) );
    else
      _os.reset( new std::ofstream( _outputFilename.c_str() ) );
      

    _solutionsToPrint.clear(); ;
    for (unsigned int i=0; i < _input.vector_variable_size("solutions"); i++)
      _solutionsToPrint.push_back(_input("solutions", "",i) );

    _computeMeans        = _input("computeMeans",true);
    _computeNorms        = _input("computeNorms",true);
    _outputIterations    = _input("iterations",false);
    _outputCoefficients  = _input("coefficients",true);
    _outputWeights       = _input("weights",true);
    _outputPoints        = _input("points",true);
    _outputIndexSet      = _input("index_set",true);

    _generateSamples     = _input("generateSamples",false);
    _sampleFile          = _input("sampleFile","./sampleFile");
    _nSamples            = _input("nSamples",10000);

    // reset prefix to root
    _input.set_prefix("") ;
  }

/********************************************//**
 * \brief build physics model from default classes.
 * Can be overidden to include new physicsModel classes.
 * 
 ***********************************************/
  std::shared_ptr<PhysicsModel<T_S,T_P> > Driver::_initPhysics( GetPot& input)
  {
    //TODO  deal with parallel issue
    //what?
    std::shared_ptr< PhysicsModel<T_S,T_P> > physics ;
    std::string physicsName = input("physics/type","catenary");
    _refinePhysics = input("physics/refine",false);
    _uniformRefine = input("physics/uniformRefine",true);

    if(AGNOS_DEBUG)
      std::cout << "_initPhysics() rank: " << _physicsComm.rank() << std::endl;

    if ( physicsName == "viscousBurgers" )
    {
      input.set_prefix("physics/viscousBurgers/");
      physics =
          std::shared_ptr<AGNOS::PhysicsViscousBurgers<T_S,T_P> >(
            new AGNOS::PhysicsViscousBurgers<T_S,T_P>(
              _physicsComm, input )
            ) ;
      input.set_prefix("");
    }
#ifdef AGNOS_ENABLE_CHANNELFLOW
    else if ( physicsName == "channelFlow" )
    {
      input.set_prefix("physics/");
      physics =
          std::shared_ptr<AGNOS::PhysicsChannelFlow<T_S,T_P> >(
            new AGNOS::PhysicsChannelFlow<T_S,T_P>( _physicsComm, input )
            ) ;
      input.set_prefix("");
    }
#endif // AGNOS_ENABLE_CHANNEL_FLOW
#ifdef AGNOS_ENABLE_GRINS
    else if ( physicsName == "grins" )
    {
      input.set_prefix("physics/");
      physics =
          std::shared_ptr<AGNOS::PhysicsGrins<T_S,T_P> >(
            new AGNOS::PhysicsGrins<T_S,T_P>( _physicsComm, input )
            ) ;
      input.set_prefix("");
    }
#endif // AGNOS_ENABLE_GRINS
    else if ( physicsName == "catenaryLibmesh" )
    {
      input.set_prefix("physics/catenaryLibmesh/");
      physics =
          std::shared_ptr<AGNOS::PhysicsCatenaryLibmesh<T_S,T_P> >(
            new AGNOS::PhysicsCatenaryLibmesh<T_S,T_P>( _physicsComm, input )
            ) ;
      input.set_prefix("");
    }
    else if ( physicsName == "catenary" )
    {
      input.set_prefix("physics/catenary/");
      physics =
          std::shared_ptr<AGNOS::PhysicsCatenary<T_S,T_P> >(
            new AGNOS::PhysicsCatenary<T_S,T_P>(
              _physicsComm, input )
            ) ;
      input.set_prefix("");
    }
    else if ( physicsName == "diffusion" )
    {
      input.set_prefix("physics/diffusion/");
      physics =
          std::shared_ptr<AGNOS::PhysicsDiffusion<T_S,T_P> >(
            new AGNOS::PhysicsDiffusion<T_S,T_P>(
              _physicsComm, input )
            ) ;
      input.set_prefix("");
    }
    else
    {
      std::cerr << " ERROR: unrecognized physics type: " << physicsName
        << "\n"
        << "please choose an appropriate physics type "
        << "(or implememnt it yourself) \n"
        << std::endl;
      exit(1);
    }

    //if exactQoi is present set the flag
    this->_exactQoiExists = ( 
        physics->getAvailableSolutions().find("exactQoi")
        != physics->getAvailableSolutions().end() 
        ) ;

    return physics;
  }

/********************************************//**
 * \brief build primary surrogateModel
 * Settings and options provided through libMesh input file. 
 *
 * handles the initialziation of each surrogate model listed in input file
 * 
 ***********************************************/
  std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > >
    Driver::_initSurrogate( GetPot& input, 
        std::vector< std::shared_ptr<AGNOS::Parameter> > parameters,
        std::shared_ptr<PhysicsModel<T_S,T_P> > physics)
  {
    input.set_prefix("surrogateModels/");
    unsigned int nSurrogates = input.vector_variable_size("modelNames");
    std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates;
    surrogates.reserve(nSurrogates);

    /* _surrogateNames.resize( input.vector_variable_size("modelNames") ); */
    for(unsigned int n=0; n < nSurrogates; n++)
    {
      std::string modelName = input("modelNames","", n);
      _surrogateNames.insert( std::pair<std::string,unsigned int>(
            modelName, n )
          );

      // get section name and set prefix
      std::string sectionName = "surrogateModels/";
        sectionName += modelName;
        sectionName += "/";
      input.set_prefix(sectionName.data()) ;


      /** Get type of surrogate model */
      std::string surrType  = input("type","PseudoSpectralTensorProduct");
      int surrogateType;

      if (surrType == "PseudoSpectralTensorProduct")
        surrogateType =  PSEUDO_SPECTRAL_TENSOR_PRODUCT ;
      else if (surrType == "PseudoSpectralSparseGrid")
        surrogateType =  PSEUDO_SPECTRAL_SPARSE_GRID ;
      else if (surrType == "PseudoSpectralMonteCarlo")
        surrogateType =  PSEUDO_SPECTRAL_MONTE_CARLO ;
      else if (surrType == "EvaluatorPseudoSpectral")
        surrogateType =  EVALUATOR_PSEUDO_SPECTRAL ;
      else if (surrType == "Collocation")
        surrogateType =  COLLOCATION ;
      else
      {
        std::cerr << " ERROR: unrecognized SurrogateModelType " 
          << surrType << std::endl;
        exit(1);
      }

      /** solutions we want the surrogate to approximate */
      std::set<std::string> computeSolutions ;
      for(unsigned int n=0; n < input.vector_variable_size("computeSolutions"); n++)
        computeSolutions.insert( input("computeSolutions","", n) );

      //if exactQoi is present set the flag
      if (n==0) // primary (or at least the first) surrogate
        this->_exactQoiExists *= ( computeSolutions.find("exactQoi") != computeSolutions.end() );


      /** Determine which type of surrogate we have, primary or secondary */
      std::string primarySurrogateName = input("primarySurrogate","");

      /** solutions we will evaluate from a primary surrogate */
      std::set<std::string> evaluateSolutions ;
      for(unsigned int n=0; n < input.vector_variable_size("evaluateSolutions"); n++)
        evaluateSolutions.insert( input("evaluateSolutions","", n) );


      /** Get length of input vectors */
      int orderDim = input.vector_variable_size("order");
      std::vector<unsigned int> order;

      if (orderDim == 1)
        for (unsigned int i=0; i < _paramDim; i++)
          order.push_back( input("order", 0) ) ;
      else
        for (unsigned int i=0; i < _paramDim; i++)
          order.push_back( input("order", 0, i) ) ;

      /** Get increaseOrder */
      /** Get length of input vectors */
      int incOrderDim = input.vector_variable_size("increaseOrder");
      std::vector<unsigned int> increaseOrder;;

      if (incOrderDim == 1)
        for (unsigned int i=0; i < _paramDim; i++)
          increaseOrder.push_back( input("increaseOrder", 0) ) ;
      else
        for (unsigned int i=0; i < _paramDim; i++)
          increaseOrder.push_back( input("increaseOrder", 0, i) ) ;
      
      /** Get multiplyOrder */
      unsigned int multiplyOrder = input("multiplyOrder",1);



      // type specific setup
      switch( surrogateType )
      {
        case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
          {
            /** primary surrogate  */
            if ( primarySurrogateName.empty() )
            {
              surrogates.push_back(
                  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
                    new PseudoSpectralTensorProduct<T_S,T_P>( 
                      _comm, 
                      physics,
                      parameters, 
                      order,
                      computeSolutions )
                    ) 
                  );
            }
            /** secondary surrogate  */
            else 
            {
              std::shared_ptr<AGNOS::SurrogateModelBase<T_S,T_P> >
                primarySurrogate =
                surrogates[_surrogateNames[primarySurrogateName]] ;
              surrogates.push_back(
                  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
                    new PseudoSpectralTensorProduct<T_S,T_P>( 
                      primarySurrogate,
                      increaseOrder,
                      multiplyOrder,
                      evaluateSolutions,
                      computeSolutions
                      )
                    ) 
                  );
            }


            break;
          }
        case(PSEUDO_SPECTRAL_SPARSE_GRID):
          {
            std::cerr 
              << " this SurrogateModelType is not yet implemented\n" ;
            exit(1);
            break;
          }
        case(PSEUDO_SPECTRAL_MONTE_CARLO):
          {
            std::cerr 
              << " this SurrogateModelType is not yet implemented\n" ;
            exit(1);
            break;
          }
        case(COLLOCATION):
          {
            std::cerr 
              << " this SurrogateModelType is not yet implemented\n" ;
            exit(1);
            break;
          }
      }
    
      if(AGNOS_DEBUG)
        std::cout << "sectionName: " << sectionName << std::endl;
      input.set_prefix("surrogateModels/") ;
    }

    input.set_prefix("");

    return surrogates;
  }

/********************************************//**
 * \brief build a simulation driver using a restart file
 * 
 ***********************************************/
 void  Driver::buildFromRestart( )
  {
    // initialize as usual
    // read in restart data
  }

/********************************************//**
 * \brief build a simulation driver using a restart file
 * 
 ***********************************************/
 void  Driver::buildEvaluator( H5IO& h5io )
  {
    _activeElems.clear();
    // read in active elements
    h5io.readSimulation( _activeElems, _comm, _physicsComm );
  }

/********************************************//**
 * \brief An evaluator routine to compute the global surrogate model
 * 
 * Needed for external evaluators, since the surrogate model may be constructed
 * of SurrogateEvaluators on multiple elements. 
 *
 * Defaults to evaluating surrogate 0 on each element and returning all
 * solutions.
 ***********************************************/
  T_P Driver::evaluate( std::string solutionName, T_S& parameterValues )
  {
    if(AGNOS_DEBUG)
    {
      std::cout 
        << "------------- Evaluating global surrogate model ----- "
        << std::endl;
      std::cout << "nActiveElems = " << _activeElems.size() << std::endl;
    }

    T_P value;

    // loop through active elements 
    std::list<AGNOS::Element<T_S,T_P> >::iterator elit ;
    bool covered = false;
    for (elit=_activeElems.begin(); elit!=_activeElems.end(); ++elit)
      if ( elit->covers(parameterValues) )
      {
        covered = true;
        value = 
          elit->surrogates()[0]->evaluate( solutionName, parameterValues );
      }

    if (!covered)
    {
      std::cerr << "***** parameter value not in support of surrogate model \n";
      std::cerr << "p: " ;
      for(unsigned int p=0;p<parameterValues.size();p++)
        std::cerr << parameterValues(p) << " " ;
      std::cerr << std::endl;
      std::cerr << " throwing error  ***********" << std::endl;
      throw -1;
    }


    if(AGNOS_DEBUG)
      std::cout << "leaving Driver::evaluate" << std::endl;
    return value;
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void Driver::run( )
  {

    std::ofstream errorOut;
    errorOut.open("error.txt", std::ios::trunc);
    errorOut << "#" << " " ;
    errorOut << "N_coeffs" << " " ;
    errorOut << "N_phy_dof" << " " ;
    if (_exactQoiExists)
      errorOut << "E_exact" << " "  ;
    errorOut << "E_total" << " "  ;
    errorOut << "E_physics" << " " ;
    errorOut << "E_surrogate" << std::endl;

    // First iteration
    // TODO can we combine this into iter loop?
    std::cout 
      << "----------------- ITER " << 1 << " -----------------" 
      << std::endl;

    while (!_elemsToUpdate.empty())
    {
      const AGNOS::Element<T_S,T_P>& elem = _elemsToUpdate.front();
      // build initial approximation
      for(unsigned int i=0;i<elem.surrogates().size(); i++)
      {
        std::cout << "pre surrogate build " << i << std::endl;
        elem.surrogates()[i]->build();
        std::cout << "post surrogate build " << i << std::endl;
      }

      // add element to active list
      _activeElems.push_front(elem) ;

      // remove element from update list
      _elemsToUpdate.pop();
    }


    std::list<AGNOS::Element<T_S,T_P> >::iterator elit =
      _activeElems.begin();

    // initialize dof counters
    double globalNCounter = 0;
    double globalSolSizeCounter = 0;
    
    // initialize global errors
    double globalExactError = 0;
    double globalTotalErrorEstimate = 0;
    double globalPhysicsErrorEstimate = 0;
    double globalSurrogateErrorEstimate = 0;
    double maxElementError = 0;

    // if we are using adaptiveDriver then print out/save individual errors
    if (_adaptiveDriver)
    {
      // loop through all active elements
      for (; elit!=_activeElems.end(); elit++)
      {

        globalNCounter        += elit->surrogates()[0]->getTotalNCoeff();
        globalSolSizeCounter  += elit->surrogates()[0]->getPhysics()->getNDofs();



        computeErrors( 
            *elit, 
            globalPhysicsErrorEstimate,
            globalTotalErrorEstimate,
            globalSurrogateErrorEstimate,
            maxElementError);

        if(_exactQoiExists)
        {
          // square current sums so we can add additional elements as sum of
          // squares
          globalExactError *= globalExactError;
          // get exact error contrib
          double localExactError = elit->surrogates()[0]->evaluateError("exactQoi");
          globalExactError += elit->weight() * std::pow(localExactError,2.) ;
          globalExactError  = std::sqrt( globalExactError );
        }


      } // end for active elements

    } // end if adaptiveDriver

    // broadcast errors to all procs
    if(_exactQoiExists)
      _physicsComm.broadcast(globalExactError);
    _physicsComm.broadcast(globalTotalErrorEstimate);
    _physicsComm.broadcast(globalSurrogateErrorEstimate);
    _physicsComm.broadcast(globalPhysicsErrorEstimate);

    std::cout << "     N_ElEM:  "  << _activeElems.size() << std::endl;
    std::cout << "   N_COEFFS:  "  << globalNCounter << std::endl;
    std::cout << " N_PHY_DOFS:  "  << globalSolSizeCounter << std::endl;
    if (_exactQoiExists)
      std::cout << "GLOBAL:  totalExactError         = "  << globalExactError << std::endl;
    std::cout << "GLOBAL:  totalErrorEstimate      = "  << globalTotalErrorEstimate << std::endl;
    std::cout << "GLOBAL:  physicsErrorEstimate    = "  << globalPhysicsErrorEstimate << std::endl;
    std::cout << "GLOBAL:  surrogateErrorEstimate  = "  << globalSurrogateErrorEstimate << std::endl;

    errorOut << globalNCounter << " " ;
    errorOut << globalSolSizeCounter << " " ;
    if (_exactQoiExists)
      errorOut << globalExactError << " "  ;
    errorOut << globalTotalErrorEstimate << " "  ;
    errorOut << globalPhysicsErrorEstimate << " " ;
    errorOut << globalSurrogateErrorEstimate << std::endl;



    // now perform the rest of the iterations
    for (unsigned int iter=2; iter <= _maxIter; iter++)
    {

    maxElementError = 0;
    std::cout 
      << "----------------- ITER " << iter << " -----------------" 
      << std::endl;

    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
    if(AGNOS_DEBUG)
      std::cout << "DEBUG: iter-" << iter << " rank-" << globalRank << std::endl;

      
      // need a way to keep track of whether physics has been refined or not
      std::set< std::shared_ptr<PhysicsModel<T_S,T_P> > > markedPhysics;
      
      // loop through active elements , if error is above threshold mark
      // unique physics objects
      for (elit=_activeElems.begin(); elit!=_activeElems.end(); ++elit)
      {
        if ( elit->_totalErrorEstimate >= _refinePercentage * maxElementError )
        {
          // Determine which space to refine
          if ( ( _refinePhysics 
              && 
              ( //(!_adaptiveDriver)
                //||
                (globalSurrogateErrorEstimate <= globalPhysicsErrorEstimate)  
              )
              ) || _simultRefine || _forcePhysicsRefine
             ) 
          {
            /* if(AGNOS_DEBUG) */
              std::cout << "DEBUG: refine physics rank-" << globalRank << std::endl;
          
            // mark my physics to be refined if it hasn't been
            if ( !markedPhysics.count( elit->physics() ) ) 
            {
              if(AGNOS_DEBUG)
                std::cout << "DEBUG: unique physics" << std::endl;
              markedPhysics.insert( elit->physics() ) ;
            }

            // TODO this may not be enough, we may need to reistantiate the
            // surrogate model, or at least set physics
            
            // add to update queue
            _elemsToUpdate.push(*elit);

          } // end of if refine physics

          //otherwise refine surrogate
          if ( ( _refineSurrogate 
              &&
              ( //(!_adaptiveDriver)
                //||
                (globalSurrogateErrorEstimate >= globalPhysicsErrorEstimate)
              )
              ) || _simultRefine || _forceSurrogateRefine
              )
          {

            /* if(AGNOS_DEBUG) */
              std::cout << "DEBUG: refine surrogate rank-" 
                << globalRank << std::endl;
             

            // test if h refinement is practical
            double childMaxSurrogateError = 0.;
            std::vector< Element<T_S,T_P> > children ;
            if ( _hRefine ) // TODO conditions
            {
              /* if(AGNOS_DEBUG) */
                std::cout << "DEBUG: h refineing element rank-" 
                  << globalRank << std::endl;

              children = elit->split() ;

              double childPhysicsError = 0.;
              double childTotalError = 0.;
              double childSurrogateError = 0.;
              double childMaxError = 0.;
              double childErrorSum = 0.;
              for (unsigned int c=0; c<children.size(); c++)
              {
                computeErrors(
                    children[c],
                    childPhysicsError,
                    childTotalError,
                    childSurrogateError,
                    childMaxError
                    );

                // add contribution to sum
                childErrorSum += childSurrogateError ;
                // if this child has larger error save as max
                if ( childSurrogateError > childMaxSurrogateError)
                  childMaxSurrogateError = childSurrogateError ;
              }

              // make sure contributions roughly equal element error
              if( std::abs(childErrorSum - elit->_totalErrorEstimate) <= 1e-9 )
              {
                std::cerr << std::endl;
                std::cerr 
                  << " ERROR: children error contributions don't add up to \n"
                  << "        total element error" 
                  << std::endl;
                std::cerr << std::endl;
                exit(1);
              } // end if contributions match
            } /// end if hRefine


            // perform h refinement
            // TODO make percentage an option??
            if ( _hRefine 
                && 
                (childMaxSurrogateError>= (0.25 * elit->_surrogateErrorEstimate) ) 
                )
            {
              
              // cosntruct new surrogates for each child element
              // TODO: optionally keep old surrogate for elems with non-dominant
              // error contributions
              for (unsigned int c=0; c<children.size(); c++)
              {
                std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > 
                  childSurrogates = _initSurrogate( 
                      _input, children[c].parameters(), children[c].physics() ) ;
                children[c].setSurrogates( childSurrogates ); 

                // save element to update queue
                _activeElems.push_front(children[c]) ;
                _elemsToUpdate.push(children[c]) ;
              }
              // remove old element from active too
              elit = _activeElems.erase(elit);
              elit--;

            }
            // perform p refinement
            else if( _pRefine )
            {
              if(AGNOS_DEBUG)
                std::cout << "DEBUG: p refineing element rank-" 
                  << globalRank << std::endl;
              std::vector<unsigned int> increase(_paramDim,0) ;
              if (!_anisotropic)
              {
              if(AGNOS_DEBUG)
                std::cout << "DEBUG: uniform p refineing element rank-" 
                  << globalRank << std::endl;
                for (unsigned int i=0; i< increase.size(); i++)
                  increase[i] += _pIncrement[i];
              }
              // if anisotropic refinement is being used we need to determine
              // which direction to refine in
              else
              {
                /* if(AGNOS_DEBUG) */
                  std::cout << "DEBUG: anisotropic p refineing element rank-" 
                    << globalRank << std::endl;
                // TODO we could also give option for other types of anisotropic
                // refinement
                
                // if adaptiveDriver, use error Surrogate
                if (_adaptiveDriver)
                {
                  // get error coefficients
                  std::map< std::string, LocalMatrix> errorCoeffs =
                    elit->surrogates()[1]->getCoefficients() ;
    

                  // copy index sets (these are already sorted by construction)
                  std::vector< std::vector<unsigned int> > sortedN 
                    = elit->surrogates()[0]->indexSet() ;
                  std::vector< std::vector<unsigned int> > sortedM 
                    = elit->surrogates()[1]->indexSet() ;

                  /* if(AGNOS_DEBUG) */
                  {
                    std::cout << " N: " << std::endl;
                    for (unsigned int i=0; i<sortedN.size(); i++)
                    {
                      for (unsigned int j=0; j<sortedN[i].size(); j++)
                        std::cout << sortedN[i][j] << " " ;
                      std::cout << std::endl;
                          
                    }
                    std::cout << " sorted M: " << std::endl;
                    for (unsigned int i=0; i<sortedM.size(); i++)
                    {
                      for (unsigned int j=0; j<sortedM[i].size(); j++)
                        std::cout << sortedM[i][j] << " " ;
                      std::cout << std::endl;
                          
                    }
                  }


                  // determine largest coeff of error surrogate
                  unsigned int i, unique, maxIndex = 0;
                  double maxCoefficient = 0; 
                  std::vector<std::vector<unsigned int> >::iterator sit = sortedM.begin(); 

                  if (!errorCoeffs.empty())
                  {
                    for (i=0,unique=0;sit!= sortedM.end(); sit++, i++)
                    {
                      if ( *sit != sortedN[i-unique] )
                      {
                        unique++;

                        if ( std::abs(errorCoeffs["errorEstimate"](i,0)) >= maxCoefficient )
                        {
                          maxCoefficient = std::abs(errorCoeffs["errorEstimate"](i,0)) ;
                          maxIndex = i;
                        }

                        /* if(AGNOS_DEBUG) */
                        {
                          std::cout << " found one: index:" << i << "    " ;
                          for (unsigned int j=0; j<sit->size(); j++)
                            std::cout << (*sit)[j] << " " ;
                          std::cout << "   coeff value = " << errorCoeffs["errorEstimate"](i,0) ;
                          std::cout << std::endl;
                        } // if AGNOS_DEBUG
                      } // if unique index
                    } // for each index M
                  } // if errorCoeffs not empty

                  this->_comm.broadcast(maxIndex,0);
                  this->_physicsComm.broadcast(maxIndex,0);

                  if(AGNOS_DEBUG)
                    std::cout << "   maxIndex = " << maxIndex << std::endl;
                  
                  // increase only that dir
                  for (unsigned int j=0; j<sortedM[maxIndex].size();j++)
                    if ( sortedM[maxIndex][j] >
                        elit->surrogates()[0]->getExpansionOrder()[j] )
                      increase[j] += _pIncrement[j];

                  if(AGNOS_DEBUG)
                  {
                    std::cout << "increase = " ;
                    for (unsigned int i=0; i<increase.size(); i++)
                      std::cout << increase[i] << " " ;
                    std::cout << std::endl;
                  }


                  
                } // end if adaptiveDriver

              } // end if anisotropic


              // refine all surrogates. We don't need to worry about the
              // anisotropic causing problems with secondary surrogates,
              // refinement based on primary takes precedence. 
              for(unsigned int i=0;i<elit->surrogates().size(); i++)
                elit->surrogates()[i]->refine( increase ) ;
              
              // add to update queue
              _elemsToUpdate.push(*elit);

            }
            else
            {
              std::cout << std::endl;
              std::cerr << 
                " ERROR: surrogate marked for refinement but both hRefine ";
              std::cerr 
                << "        and pRefine are set to false. "
                << std::endl;
              std::cout << std::endl;
              exit(1);
            }
            
          } // end of if refine surrogate

        } // end of active element and if above threshold loop
        else
        {
          if(AGNOS_DEBUG)
          {
            std::cout << "DEBUG: no refine rank-" << globalRank << std::endl;
            std::cout << "      totalErrorEstimate=" << elit->_totalErrorEstimate ;
            std::cout << "      maxElementError=" << maxElementError;
            std::cout << std::endl;
          }
        }

        if ( elit == _activeElems.end() )
          break;
      } // loop over active elements
      std::cout << "nActiveElems = " << _activeElems.size() << std::endl;

      /* if(AGNOS_DEBUG) */
        std::cout << "DEBUG: pre update elements rank-" << globalRank << std::endl;

      // refine all physics that were marked
      std::set< std::shared_ptr<PhysicsModel<T_S,T_P> > >::iterator physIt =
        markedPhysics.begin();
      for(; physIt != markedPhysics.end(); physIt++)
      {
        if(_uniformRefine)
          (*physIt)->refine()  ;
        else
        {
          // compute means of errorIndicators
          // if they aren't present in the surrogate model, assertion will be
          // throwin inside computeMeans
          std::map<std::string,T_P> indicatorMeans;
          std::vector<std::string> indName(1,"errorIndicators");
          computeMeans( indName, _activeElems, indicatorMeans );
          (*physIt)->refine( indicatorMeans["errorIndicators"] );
        }
      }



      while (!_elemsToUpdate.empty())
      {
        const AGNOS::Element<T_S,T_P>& elem = _elemsToUpdate.front();

        // rebuild surrogates for elements that need updated
        for(unsigned int i=0;i<elem.surrogates().size(); i++)
        {
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: pre surrogate build " << i << " rank-" << globalRank << std::endl;
          elem.surrogates()[i]->build();
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: post surrogate build " << i << " rank-" << globalRank << std::endl;
        }

        // remove element from update list
        _elemsToUpdate.pop();
      }

      if(AGNOS_DEBUG)
        std::cout << "DEBUG: post update elements rank-" << globalRank << std::endl;
        
      // if we are using adaptiveDriver then print out/save individual errors
      if (_adaptiveDriver)
      {
        // reset global values to 0
        globalNCounter = 0.;
        globalSolSizeCounter = 0.;
        globalExactError = 0;
        globalPhysicsErrorEstimate = 0. ;
        globalSurrogateErrorEstimate = 0.;
        globalTotalErrorEstimate = 0.;

        // loop through all active elements
        std::list<AGNOS::Element<T_S,T_P> >::iterator elit =
          _activeElems.begin();
        for (; elit!=_activeElems.end(); ++elit)
        {
          globalNCounter += elit->surrogates()[0]->getTotalNCoeff();
          globalSolSizeCounter  += elit->surrogates()[0]->getPhysics()->getNDofs();

          computeErrors( 
              *elit, 
              globalPhysicsErrorEstimate,
              globalTotalErrorEstimate,
              globalSurrogateErrorEstimate,
              maxElementError);

          if(_exactQoiExists)
          {
            // square current sums so we can add additional elements as sum of
            // squares
            globalExactError *= globalExactError;
            // get exact error contrib
            double localExactError = elit->surrogates()[0]->evaluateError("exactQoi");
            globalExactError += elit->weight() * std::pow(localExactError,2.) ;
            globalExactError  = std::sqrt( globalExactError );
          }

        } // end for active elements

      } // end if adaptiveDriver
      

      // broadcast errors to all procs
      if(_exactQoiExists) 
        _physicsComm.broadcast(globalExactError);
      _physicsComm.broadcast(globalTotalErrorEstimate);
      _physicsComm.broadcast(globalSurrogateErrorEstimate);
      _physicsComm.broadcast(globalPhysicsErrorEstimate);
      
      std::cout << "     N_ElEM:  "  << _activeElems.size() << std::endl;
      std::cout << "   N_COEFFS:  "  << globalNCounter << std::endl;
      std::cout << " N_PHY_DOFS:  "  << globalSolSizeCounter << std::endl;
      if (_exactQoiExists)
        std::cout << "GLOBAL:  totalExactError         = "  << globalExactError << std::endl;
      std::cout << "GLOBAL:  totalErrorEstimate      = "  << globalTotalErrorEstimate << std::endl;
      std::cout << "GLOBAL:  physicsErrorEstimate    = "  << globalPhysicsErrorEstimate << std::endl;
      std::cout << "GLOBAL:  surrogateErrorEstimate  = "  << globalSurrogateErrorEstimate << std::endl;

      errorOut << globalNCounter << " " ;
      errorOut << globalSolSizeCounter << " " ;
      if (_exactQoiExists)
        errorOut << globalExactError << " "  ;
      errorOut << globalTotalErrorEstimate << " "  ;
      errorOut << globalPhysicsErrorEstimate << " " ;
      errorOut << globalSurrogateErrorEstimate << std::endl;


    } // end for nIter

    // close error file
    errorOut.close();


    postProcess();
    

    // save driver info to hdf5 file for restarts and evaluators
    std::cout << "    writing final simulation to file " << std::endl;
    _h5io->writeSimulation( _activeElems );


    return;
  }
  
/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::postProcess(   ) 
  {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout 
      << "----------------- POST PROCESSING ---------------------" ;
    std::cout << std::endl;

    // print out settings
    std::cout << "    printing settings " << std::endl;
    printSettings();
    std::cout << "    printing final solution data " << std::endl;
    printSolution(_maxIter);

    // compute means of final surrogate
    if(_computeMeans)
    {
      std::cout << "    computing means for final surrogate " << std::endl;
      std::ostream& out = *(this->_os) ;
      out << std::endl;
      out << "#====================================================" 
        << std::endl;
      out << "#--------- MEANS ------------------\n";

      std::map<std::string,T_P> globalMeans;
      computeMeans( _solutionsToPrint, _activeElems, globalMeans );

      for(unsigned int i=0; i<_solutionsToPrint.size();i++)
      {
        out << "solution: " << _solutionsToPrint[i] << std::endl;
        globalMeans[_solutionsToPrint[i]].print(out) ; 
      }
    }

    // compute norms of final surrogate
    if(_computeNorms)
    {
      std::cout << "    computing means for final surrogate " << std::endl;
      std::ostream& out = *(this->_os) ;
      out << std::endl;
      out << "#====================================================" 
        << std::endl;
      out << "#--------- NORMS ------------------\n";

      std::map<std::string,T_P> globalNorms;
      computeNorms( _solutionsToPrint, _activeElems, globalNorms );

      for(unsigned int i=0; i<_solutionsToPrint.size();i++)
      {
        out << "solution: " << _solutionsToPrint[i] << std::endl;
        globalNorms[_solutionsToPrint[i]].print(out) ; 
      }
    }
    
    // TODO generate samples and estimate stats 
    // sample surrogate
    std::vector<T_P> sampleVec;
    if(_generateSamples)
    {
      std::cout << "    generating samples from final surrogate " << std::endl;
      // open output file
      std::ofstream sampleOut;
      sampleOut.open(_sampleFile, std::ios::trunc);
      sampleOut << std::setprecision(5) << std::scientific ;
      
      //clear out any old samples
      sampleVec.clear();

      // make sure solutionName is available
      /* agnos_assert( (_solutionNames.count( solutionName )) ) ; */

      // variables needed for gsl randum number generator
      std::vector<const gsl_rng_type *> T;
      std::vector<gsl_rng *> r;
      gsl_rng_env_setup();

      // loop through parameters and set up rng
      for(unsigned int p=0; p<_paramDim; p++)
      {
        // initialize rng for this parameter type
        T.push_back( gsl_rng_default );
        r.push_back( gsl_rng_alloc(T.back()) );
      } // end 

      // sample generation loop
      for(unsigned int i=0; i<_nSamples; i++)
      {
        std::vector<double> s;
        for(unsigned int p=0; p<_paramDim; p++)
        {
          double scaledSample;
          // check type 
          // UNIFROM
          switch( ParameterType(_globalParameters[p]->type()) )
          {
            case UNIFORM:
              scaledSample = _globalParameters[p]->min() 
                + (_globalParameters[p]->max() - _globalParameters[p]->min() ) 
                * gsl_rng_uniform(r[p]) ;
              s.push_back( scaledSample  );
              break;
            case CONSTANT:
              s.push_back( _globalParameters[p]->min() ) ;
              break;
          }
        }

        // store as template vector type
        T_S paramValues(s);

        // evaluate model
        sampleVec.push_back( evaluate( "qoi", paramValues) ) ;


      } // end loop over samples

      // free up memory
      for(unsigned int p=0; p<_paramDim; p++)
        gsl_rng_free(r[p]);








      /* std::list<AGNOS::Element<T_S,T_P> >::iterator elit = */
      /*   _activeElems.begin(); */
      /* for (; elit!=_activeElems.end(); elit++) */
      /* { */
      /*   // at this point only sample from primary surrogate */
      /*   elit->surrogates()[0]->sample( "qoi", _nSamples, sampleVec  ); */

        for(unsigned int s=0; s<_nSamples; s++)
          sampleOut << sampleVec[s](0) << std::endl;
      /* } */

      sampleOut.close();
    }

    return;
  }



/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printDriverSettings(   ) 
  {
    std::ostream& out = *(this->_os) ;

    out << std::endl;
    out << "#====================================================" <<
      std::endl;
    out << "# Driver settings: " << std::endl;
    out << "#    maxIter = " << _maxIter << std::endl;

    out << std::endl;
    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printParameterSettings(  ) 
  {
    std::ostream& out = *(this->_os) ;
    out << std::endl;
    out << "#====================================================" <<
      std::endl;
    out << "# Parameter settings: " << std::endl;
    out << "#     dimension = " << _paramDim << std::endl;
    out << "#     nInitialHRefinements = " << _nInitialHRefinements << std::endl;
    out << "#     nElems = " <<
      std::distance(_activeElems.begin(),_activeElems.end()) << std::endl;


    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSurrogateSettings(  ) 
  {
    return;
  }


/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSolutionData(  ) 
  {
    std::ostream& out = *(this->_os) ;

    // Per surrogate options
    std::list<AGNOS::Element<T_S,T_P> >::iterator elit =
      _activeElems.begin();
    for (; elit!=_activeElems.end(); elit++)
    {
      if (this->_comm.rank() == 0)
      {
        out << std::endl;
        out << "#====================================================" 
          << std::endl;
        out << "#--------- ELEM ------------------\n";
        out << "#     mins = " ;
        for (unsigned int i=0; i < _paramDim; i++)
          out << elit->parameters()[i]->min() << " " ;
        out << std::endl;

        out << "#     maxs = " ;
        for (unsigned int i=0; i < _paramDim; i++)
          out << elit->parameters()[i]->max() << " " ;
        out << std::endl;
        out << std::endl;
      }

      std::map<std::string,unsigned int>::iterator sid =
        _surrogateNames.begin();
      for(; sid !=_surrogateNames.end(); ++sid)
      {
        if (this->_comm.rank() == 0)
        {
          out << "#----------------------------------------------------" <<
            std::endl;
          out << "#     " << sid->first << ": " << std::endl;
          out << "#     order = " ;
          std::vector<unsigned int> order = elit->surrogates()[sid->second]->getExpansionOrder();
          for(unsigned int i=0; i < _paramDim; i++)
            out << order[i] << " " ;
          out << std::endl;
          out << "#     coefficients: " << std::endl;
        }

        if (_outputCoefficients)
          elit->surrogates()[sid->second]->printCoefficients( _solutionsToPrint, out );
        if (_outputWeights)
          elit->surrogates()[sid->second]->printIntegrationWeights( out );
        if (_outputPoints)
          elit->surrogates()[sid->second]->printIntegrationPoints( out );
        if (_outputIndexSet)
          elit->surrogates()[sid->second]->printIndexSet( out );
      }

      if (this->_comm.rank() == 0)
      {
        out << "#----------------------------------------------------" <<
          std::endl;
      }
    }

    return;
  }



/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSettings( ) 
  {

      printDriverSettings( ) ;
      printParameterSettings(  ) ;

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSolution( unsigned int iteration ) 
  {
      std::ostream& out = *(this->_os) ;
      out << std::endl;
      out << "#====================================================\n" 
         << "#      Solution data for ITER " << iteration << std::endl;
      printSolutionData( ) ;

    return;
  }


/********************************************//**
 * \brief compute element errors
 ***********************************************/
  void Driver::computeErrors( 
      AGNOS::Element<T_S,T_P>&  elem,
      double& globalPhysics,
      double& globalTotal,
      double& globalSurrogate,
      double& maxElementError
      ) 
  {
    // square current sums so we can add additional elements as sum of squares
    globalPhysics *= globalPhysics;
    globalTotal *= globalTotal;
    globalSurrogate *= globalSurrogate ;

    // get physics error contrib
    double physicsErrorEstimate;
    physicsErrorEstimate = (elem.surrogates()[0]->l2Norm("errorEstimate"))(0);

    /**Add to global tally  */
    globalPhysics+= elem.weight() * std::pow(physicsErrorEstimate,2.) ;


    // safe guard against there not being a secondary surrogate
    if (elem.surrogates().size() < 2)
    {
      std::cout << std::endl;
      std::cerr << 
        " ERROR: secondary 'error' surrogate has not been constructed"
        << std::endl;
      std::cout << std::endl;
      exit(1);
    }
    else
    {
      // total error contribution
      double totalErrorEstimate     = (elem.surrogates()[1]->l2Norm("errorEstimate"))(0);
      elem._totalErrorEstimate = totalErrorEstimate ;
      globalTotal += elem.weight() * std::pow(totalErrorEstimate,2.);

      // surrogate error contribution
      double surrogateErrorEstimate = elem.surrogates()[0]->l2NormDifference( 
            *(elem.surrogates()[1]), "errorEstimate");
      elem._surrogateErrorEstimate = surrogateErrorEstimate ;
      globalSurrogate += elem.weight() * std::pow(surrogateErrorEstimate,2.) ;


      // keep track of max of error
      if (elem._totalErrorEstimate >= maxElementError)
        maxElementError = elem._totalErrorEstimate ;

    } // end if errorSurrogate exists

    // take square root of current sums 
    globalPhysics = std::sqrt( globalPhysics );
    globalTotal  = std::sqrt( globalTotal );
    globalSurrogate  = std::sqrt( globalSurrogate );

  }

  
}
