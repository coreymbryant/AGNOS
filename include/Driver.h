
#ifndef DRIVER_H
#define DRIVER_H

/* #include "agnosDefines.h" */
#include "PseudoSpectralTensorProduct.h"
/* #include "PseudoSpectralMonteCarlo.h" */
/* #include "PseudoSpectralSparseGrid.h" */
#include "PhysicsViscousBurgers.h"
#include "PhysicsCatenary.h"
#include "PhysicsCatenaryLibmesh.h"
#include "PhysicsModel.h"
#include "PhysicsLibmesh.h"
#include "Element.h"




namespace AGNOS
{

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

      void run( ) ;

      void printSolutionData( std::ostream& out ) ;
      void printSurrogateSettings( std::ostream& out ) ;
      void printParameterSettings( std::ostream& out ) ;
      void printDriverSettings( std::ostream& out  ) ;
      void printSettings( ) ;
      void printSolution( unsigned int iteration=1 ) ;

    protected:
      std::shared_ptr< PhysicsModel<T_S,T_P> > _initPhysics( GetPot& input );
      std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > > 
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
      // ---------------------
      
      // ---------------------
      // ELEMENT VARIABLES
      std::forward_list<AGNOS::Element<T_S,T_P> > _activeElems;
      std::queue<AGNOS::Element<T_S,T_P> >        _elemsToUpdate;
      // ---------------------
      

      // ---------------------
      // PARAMETERS VARIABLES
      unsigned int              _paramDim;
      /** number of times to initially h refine the parameter domain before
       * doing any computation */
      unsigned int              _nInitialHRefinements ;
      // ---------------------
      
      // ---------------------
      // SURROGATE VARIABLES
      std::map<std::string, unsigned int> _surrogateNames;
      bool _refineSurrogate;
      // ---------------------

      // ---------------------
      // PHYSICS VARIABLES
      bool _refinePhysics;
      // ---------------------

      // ---------------------
      // OUTPUT VARIABLES
      // TODO control on per surrogate basis
      std::string               _outputFilename; 
      std::vector<std::string>  _solutionsToPrint ;
      bool                      _outputIterations  ;
      bool                      _outputCoefficients  ;
      bool                      _outputErrorCoefficients  ;
      bool                      _outputWeights       ;
      bool                      _outputPoints        ;
      bool                      _outputIndexSet      ;
      // ---------------------
  };


/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver( 
      const Communicator& comm,
      const Communicator& physicsComm,
      GetPot& input 
      ) :
    _comm(comm), _physicsComm(physicsComm) 
  {

    if(AGNOS_DEBUG)
      std::cout << "rank: " << _comm.rank() <<  std::endl;

    
    // DRIVER SETTINGS
    _maxIter = input("driver/maxIter",1);
    _adaptiveDriver = input("driver/adaptive",false);
    
    // ADAPTIVE SETTINGS
    
    
    // PARAMETER SETTINGS
    _paramDim = input("parameters/dimension", 1);
    _nInitialHRefinements = input("parameters/nInitialHRefinements", 1);

    std::vector< std::shared_ptr<AGNOS::Parameter> > parameters;
    parameters.reserve(_paramDim);
    for (unsigned int i=0; i < _paramDim; i++)
      parameters.push_back( std::shared_ptr<AGNOS::Parameter>(
            new AGNOS::Parameter( 
              input("parameters/types",0,i),
              input("parameters/mins",-1.0,i), 
              input("parameters/maxs", 1.0,i))
            )
          );
    

    
    // PHYSICS SETTINGS
    std::shared_ptr< PhysicsModel<T_S,T_P> > physics = _initPhysics( input );


    // SURROGATE MODEL(s) SETTINGS
    input.set_prefix("surrogateModels/") ;
    _refineSurrogate = input("refine",false);

    /* initialize surrogate model container */
    std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > > surrogates =
      _initSurrogate( input, parameters, physics ); 


    // Construct a single initial element and add it to the update queue
    AGNOS::Element<T_S,T_P> baseElement(
          _comm,
          parameters,
          surrogates,
          physics
          ) ;


    //TODO perform initial h refinements
    _elemsToUpdate.push( baseElement );

    

    // OUTPUT DATA SETTINGS
    input.set_prefix("output/") ;
    _outputFilename      = input("filename","cout");

    _solutionsToPrint.clear(); ;
    for (unsigned int i=0; i < input.vector_variable_size("solutions"); i++)
      _solutionsToPrint.push_back(input("solutions", "",i) );

    _outputIterations    = input("iterations",false);
    _outputCoefficients  = input("coefficients",true);
    _outputWeights       = input("weights",true);
    _outputPoints        = input("points",true);
    _outputIndexSet      = input("index_set",true);

  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::~Driver( )
  {
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
    std::string physicsName = input("physics/type","");
    _refinePhysics = input("refine",false);

    if(AGNOS_DEBUG)
      std::cout << "_initPhysics() rank: " << _comm.rank() << std::endl;

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
    else
    {
      std::cerr << " ERROR: unrecognized physics type: " << physicsName
        << "\n"
        << "please choose an appropriate physics type "
        << "(or implememnt it yourself) \n"
        << std::endl;
      exit(1);
    }

    return physics;
  }

/********************************************//**
 * \brief build primary surrogateModel
 * Settings and options provided through libMesh input file. 
 *
 * handles the initialziation of each surrogate model listed in input file
 * 
 ***********************************************/
  std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > >
    Driver::_initSurrogate( GetPot& input, 
        std::vector< std::shared_ptr<AGNOS::Parameter> > parameters,
        std::shared_ptr<PhysicsModel<T_S,T_P> > physics)
  {
    input.set_prefix("surrogateModels/");
    unsigned int nSurrogates = input.vector_variable_size("modelNames");
    std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > > surrogates;
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

      /** solutions we will evaluate from a primary surrogate */
      std::set<std::string> evaluateSolutions ;
      for(unsigned int n=0; n < input.vector_variable_size("evaluateSolutions"); n++)
        evaluateSolutions.insert( input("evaluateSolutions","", n) );


      /** Determine which type of surrogate we have, primary or secondary */
      std::string primarySurrogate = input("primarySurrogate","");

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
      unsigned int increaseOrder = input("increaseOrder",0);
      /** Get increaseOrder */
      unsigned int multiplyOrder = input("multiplyOrder",1);



      // type specific setup
      switch( surrogateType )
      {
        case(PSEUDO_SPECTRAL_TENSOR_PRODUCT):
          {
            /** primary surrogate  */
            if ( primarySurrogate.empty() )
            {
              surrogates.push_back(
                  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
                    new PseudoSpectralTensorProduct<T_S,T_P>( 
                      _comm, 
                      physics,
                      parameters, 
                      order  )
                    ) 
                  );
            }
            /** secondary surrogate  */
            else 
            {
              surrogates.push_back(
                  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
                    new PseudoSpectralTensorProduct<T_S,T_P>( 
                      surrogates[_surrogateNames[primarySurrogate]].get(),
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
    
      std::cout << "sectionName: " << sectionName << std::endl;
      input.set_prefix("surrogateModels/") ;
    }

    input.set_prefix("");


    return surrogates;
  }


/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void Driver::run( )
  {

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

      // build error surrogate
      /* _errorSurrogate->build(); */

      std::map< std::string, LocalMatrix > myCoeff =
        elem.surrogates()[0]->getCoefficients();

      // print out settings
      printSettings();
      printSolution(1);
      _elemsToUpdate.pop();
    }
    

    
    /* // print out first iteration if requested */
    /* if (this->_outputIterations && (_comm.rank() == 0) ) */
    /* { */
    /*   std::cout << "\n writing results to: " << this->_outputFilename */
    /*     << " (iter = " << 1 << " )" */
    /*    << std::endl; */
    /*   std::cout << std::endl; */
    /*   printSolution(1); */
    /* } */
    
    /*   // evaluate QoI */
    /*   T_S evalPoint(2); */
    /*   evalPoint(0) = 0.5; */
    /*   evalPoint(1) = 0.5; */
    /*   /1* T_S evalPoint(1); *1/ */
    /*   /1* evalPoint(0) = 1.5; *1/ */
    /*   T_P solutionVec = _surrogate->evaluate("primal", evalPoint ); */
    /*   T_P qoiValue = _physics[0]->evaluateQoi( evalPoint, solutionVec ) ; */

    /*   T_P l2normofphyerror = _surrogate->l2Norm("error"); */
    /*   T_P l2normoftotalerror = _errorSurrogate->l2Norm("error"); */
    /*   double normDiff = _surrogate->l2NormDifference( *_errorSurrogate, "error"); */

    /* if (_comm.rank() == 0) */
    /* { */

    /*   std::cout << "phyErrorNorm = " << l2normofphyerror << std::endl; */
    /*   std::cout << "totalErrorNorm = " << l2normoftotalerror << std::endl; */
    /*   std::cout << "| phyErrorNorm-totalErrorNorm | = " << */
    /*     std::abs(l2normofphyerror(0)-l2normoftotalerror(0)) << std::endl; */
    /*   std::cout << "normDiff = " << normDiff << std::endl; */

    /*   std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
    /* } */

    /* // refine approximation */
    /* for (unsigned int iter=2; iter <= this->_maxIter; iter++) */
    /* { */
    /*   if (_comm.rank() == 0) */ 
    /*     std::cout << "\n-------------  ITER " << iter << "  -------------\n " ; */
      

    /*   // if physical error dominates and we are allowed to refine physical */
    /*   // solution */
    /*   if ( */ 
    /*       ( _refinePhysics && (l2normofphyerror(0) >= normDiff) ) */ 
    /*       || */
    /*       ( _refinePhysics && !m_refineSurrogate ) */
    /*       ) */
    /*   { */
    /*     if (_comm.rank() == 0) */ 
    /*       std::cout << "    refining physical solution " << std::endl; */
        
    /*     // if using uniform refinement we won't have error indicators in the */
    /*     // surrogate model (by design) */
    /*     if ( _physicsFunctions.count("indicators") == 0 ) */
    /*       _physics[0]->refine( ); */
    /*     else */
    /*     { */
    /*       // retrieve first coefficient (i.e. the mean) of error inidcators */
    /*       T_P errorIndicators = (_surrogate->mean())["indicators"]; */
          
    /*       _physics[0]->refine( errorIndicators ); */
    /*     } */
    /*     if (_comm.rank() == 0) */ 
    /*       std::cout << "-------------------------------------\n " ; */

    /*     _surrogate->build(); */


    /*     _errorSurrogate->build(); */
    /*   } */


    /*   // if surrogate error dominates and we are allowed to refine surrogate */
    /*   // model */
    /*   if ( */ 
    /*       ( _refineSurrogate && (l2normofphyerror(0) < normDiff) ) */ 
    /*       || */
    /*       ( _refineSurrogate && !m_refinePhysical ) */ 
    /*       ) */
    /*   { */
    /*     if (_comm.rank() == 0) */ 
    /*     { */
    /*       std::cout << "    refining surrogate model " << std::endl; */
    /*       std::cout << "-------------------------------------\n " ; */
    /*     } */
    /*     _surrogate->refine( ); */
    /*     _errorSurrogate->refine( ); */
    /*   } */





    /*   if (this->_outputIterations && (_comm.rank() == 0) ) */
    /*   { */
    /*     std::cout << "\n writing results to: " << this->_outputFilename */
    /*       << " (iter = " << iter << " )" */
    /*       << std::endl; */
    /*     std::cout << std::endl; */
    /*     printSolution(iter); */
    /*   } */

    /*   solutionVec = _surrogate->evaluate("primal", evalPoint ); */
    /*   qoiValue = _physics[0]->evaluateQoi( evalPoint, solutionVec ) ; */

    /*   l2normofphyerror = _surrogate->l2Norm("error"); */
    /*   l2normoftotalerror = _errorSurrogate->l2Norm("error"); */
    /*   normDiff = _surrogate->l2NormDifference( *_errorSurrogate, "error"); */

    /*   if (_comm.rank() == 0) */
    /*   { */

    /*     std::cout << "phyErrorNorm = " << l2normofphyerror << std::endl; */
    /*     std::cout << "totalErrorNorm = " << l2normoftotalerror << std::endl; */
    /*     /1* std::cout << "| phyErrorNorm-totalErrorNorm | = " << *1/ */
    /*     /1* std::abs(l2normofphyerror(0)-l2normoftotalerror(0)) << std::endl; *1/ */
    /*     std::cout << "normDiff = " << normDiff << std::endl; */

    /*     std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
    /*   } */

      
      /* // evaluate QoI */
      /* T_S evalPoint(1); */
      /* evalPoint(0) = 1.5; */
      /* T_P solutionVec = _surrogate->evaluate("primal", evalPoint ); */
      /* T_P qoiValue = _physics->evaluateQoi( evalPoint, solutionVec ) ; */
      /* if (_comm.rank() == 0) */
      /* { */
      /*   std::cout << "\n Qoi = " << qoiValue(0) << std::endl; */
      /* } */
    /* } */
    
    /* // output whatever user asks for */
    /* if (_comm.rank() == 0) */
    /* { */
    /*   std::cout << "\n writing final results to: " << this->_outputFilename */
    /*     << std::endl; */
    /*   std::cout << std::endl; */
    /* } */
    /*   printSolution(_maxIter); */






    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printDriverSettings( std::ostream& out  ) 
  {
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
  void Driver::printParameterSettings( std::ostream& out ) 
  {
    std::forward_list<AGNOS::Element<T_S,T_P> >::iterator elit =
      _activeElems.begin();
    for (; elit!=_activeElems.end(); elit++)
    {
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Parameter settings: " << std::endl;
      out << "#     dimension = " << _paramDim << std::endl;

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

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSurrogateSettings( std::ostream& out ) 
  {

    out << std::endl;
    out << "#====================================================" <<
      std::endl;
    out << "# Surrogate settings: " << std::endl;

    std::map<std::string,unsigned int>::iterator sid = _surrogateNames.begin();
    for(; sid !=_surrogateNames.end(); sid++)
    {
      std::forward_list<AGNOS::Element<T_S,T_P> >::iterator elit =
        _activeElems.begin();
      for (; elit!=_activeElems.end(); elit++)
      {
          out << "#----------------------------------------------------" <<
          std::endl;
        out << "#     " << sid->first << ": " << std::endl;
        out << "#     order = " ;
        std::vector<unsigned int> order = elit->surrogates()[sid->second]->getExpansionOrder();
        for(unsigned int i=0; i < _paramDim; i++)
          out << order[i] << " " ;
        out << std::endl;
        out << "#----------------------------------------------------" <<
          std::endl;
      }
    }
    out << std::endl;
    return;
  }


/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSolutionData( std::ostream& out ) 
  {
    // Per surrogate options
    if (this->_comm.rank() == 0)
    {
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Surrogate Coefficients: " << std::endl;
    }

    std::map<std::string,unsigned int>::iterator sid = _surrogateNames.begin();
    for(; sid !=_surrogateNames.end(); ++sid)
    {
      if (this->_comm.rank() == 0)
      {
        out << "#----------------------------------------------------" <<
          std::endl;
        out << "#     " << sid->first << ": " << std::endl;
      }

      std::forward_list<AGNOS::Element<T_S,T_P> >::iterator elit =
        _activeElems.begin();
      for (; elit!=_activeElems.end(); elit++)
      {
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

    // set output steam
    if (_outputFilename == "cout" )
    {
      printDriverSettings( std::cout  ) ;
      printParameterSettings( std::cout ) ;
      printSurrogateSettings( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( _outputFilename.c_str() );

      printDriverSettings( out  ) ;
      printParameterSettings( out ) ;
      printSurrogateSettings( out ) ;

      out.close( );
    }

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  void Driver::printSolution( unsigned int iteration ) 
  {
    // set output steam
    if (_outputFilename == "cout" )
    {
      std::cout << std::endl;
      std::cout << "#====================================================\n" 
         << "#      Solution data for ITER " << iteration << std::endl;
      printSolutionData( std::cout ) ;
    }
    else
    {
      std::ofstream out;
      out.open( _outputFilename.c_str(), std::ofstream::out | std::ofstream::app );

      out << std::endl;
      out << "#====================================================\n" 
         << "#       Solution data for ITER " << iteration << std::endl;
      printSolutionData( out ) ;

      out.close( );
    }
    return;
  }
    
}
#endif // DRIVER_H


