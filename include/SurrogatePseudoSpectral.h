
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
          const std::vector<unsigned int>&  order
          );


      /** Secondary Constructor. 
       *  ***DO NOT MISTAKE FOR A COPY CONSTRUCTOR ***
       *  Intended use is for constructing a secondary surrogate model using the
       *  primary model as an evaluating object in the build routine. 
       *  If additional inputs are defined it will
       * construct a new surrogate increasing the order and using
       * primarySurrogate to perform evaluations in the constructions */
      SurrogatePseudoSpectral( 
          const SurrogateModel<T_S,T_P>* primarySurrogate, 
          unsigned int increaseOrder = 0,
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
      const std::vector< std::vector< unsigned int> > 
                                getIndexSet( ) const;
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



/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S,T_P>( comm, physics, parameters, order) 
    {
    }


/********************************************//*
 * \brief Secondary constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const SurrogateModel<T_S,T_P>* primarySurrogate, 
        unsigned int increaseOrder ,
        unsigned int multiplyOrder ,
        std::set<std::string> evaluateSolutions,
        std::set<std::string> computeSolutions
        )
      : SurrogateModel<T_S,T_P>(primarySurrogate, increaseOrder, multiplyOrder,
          evaluateSolutions, computeSolutions)
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::~SurrogatePseudoSpectral( )
    {
    }

/********************************************//**
 * \brief Get number of integration points
 ***********************************************/
  template<class T_S, class T_P>
    unsigned int SurrogatePseudoSpectral<T_S,T_P>::getNIntegrationPoints( )
    const
  {
    return _nIntegrationPoints;
  }

/********************************************//**
 * \brief Get integration weights
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<double> SurrogatePseudoSpectral<T_S,T_P>::getIntegrationWeights( )
    const
  {
    return _integrationWeights;
  }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<T_S> SurrogatePseudoSpectral<T_S,T_P>::getIntegrationPoints( )
    const
  {
    return _integrationPoints;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const std::vector< std::vector< unsigned int> > 
    SurrogatePseudoSpectral<T_S,T_P>::getIndexSet( ) const
    {
      return _indexSet;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    std::vector<double> SurrogatePseudoSpectral<T_S,T_P>::evaluateBasis( 
        const std::vector< std::vector< unsigned int > >& indexSet,
        T_S& parameterValue 
        ) const
    {
      unsigned int nTerms = indexSet.size() ;
      std::vector<double> basisValues( nTerms ,1.);

      for(unsigned int id=0; id < nTerms ; id++)
        for(unsigned int dir=0; dir < this->_dimension; dir++)
        {
          basisValues[id] 
            *= this->_parameters[dir]->evalBasisPoly( 
                indexSet[id][dir], parameterValue(dir) ) ;
        }

      return basisValues;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogatePseudoSpectral<T_S,T_P>::build( )
    {
      // This is separated from the routine that actually computes contribution
      // so that we can group surrogate models together later and wll that needs
      // to be defined is build routine based on computeContribution( )
      //
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::set< std::string >::iterator id;
      
      // output progress
      if( this->_comm.rank() == 0)
      {
        std::cout << std::endl;
        std::cout << " Starting build of surrogate model for: { " ;
        for (id=this->_solutionNames.begin();
            id!=this->_solutionNames.end(); id++)
        {
          if (id!=this->_solutionNames.begin())
            std::cout << " , " ;
          std::cout << *id;
        }
        std::cout << " }" << std::endl;
      }

      // ------------------------------------
      // Int points and coeff are currently separated in case we want to
      // implememnt higher order integration rule in the future
      // assign the appropriate int points to each processor
      unsigned int minPts = _nIntegrationPoints / this->_comm.size() ;
      unsigned int remPts = _nIntegrationPoints % this->_comm.size() ;
      unsigned int nPts  =  minPts + ( this->_comm.rank() < remPts ) ;
      unsigned int intPtsStart = this->_comm.rank()*(minPts) 
        + std::min( this->_comm.rank(), remPts ) ;

      _integrationIndices.clear();
      for (unsigned int i=0; i < nPts; i++)
        _integrationIndices.push_back( intPtsStart + i );

      // assign the appropriate coeff to each processor
      unsigned int minCoeff = totalNCoeff / this->_comm.size() ;
      unsigned int remCoeff = totalNCoeff % this->_comm.size() ;
      unsigned int nCoeffs  =  minCoeff + ( this->_comm.rank() < remCoeff ) ;
      unsigned int coeffStart = this->_comm.rank()*(minCoeff) 
        + std::min( this->_comm.rank(), remCoeff ) ;

      this->_coeffIndices.clear();
      for (unsigned int i=0; i < nCoeffs; i++)
        this->_coeffIndices.push_back( coeffStart + i );
      // ------------------------------------


      // ------------------------------------
      // Initialize data structures
      std::shared_ptr<DistMatrix> polyValues(new DistMatrix(this->_comm));
      polyValues->init( 
          this->_totalNCoeff,           // global dim m
          this->_nIntegrationPoints,    // global dim n
          nCoeffs,                      // local dim m
          nPts,                         // local dim n
          nPts,                         // # non-zero/row on proc
          this->_nIntegrationPoints-nPts// # non-zero/row off proc
          );
      polyValues->zero();

      if ( nCoeffs > 0)
      {
        assert( polyValues->row_start() == this->_coeffIndices.front() );
        assert( polyValues->row_stop() == this->_coeffIndices.back()+1  );
      }
      // ------------------------------------
      
      
      
      // ------------------------------------
      // Get primalSurrogate evaluations if needed
      this->_primaryEvaluations.clear();
      if ( ! this->_evalNames.empty() )
        for(unsigned int pt=0; pt < _nIntegrationPoints; pt++)
        {
          if (AGNOS_DEBUG)
            std::cout << "DEBUG: evaluating primarySurrogate at pt: " << pt 
              << std::endl;

          T_S integrationPoint = _integrationPoints[pt];
          bool saveLocal = ( 
              (pt >= this->_integrationIndices.front()) 
              && 
              (pt <= this->_integrationIndices.back()) 
              ) ;
          // if we don't save locally we will just be pushing an empty vector so
          // memory shouldn't be an issue
          std::map<std::string,T_P> result = this->_evalSurrogate->evaluate(
              this->_evalNames, integrationPoint, saveLocal  ) ;

          if (saveLocal)
            this->_primaryEvaluations.push_back( result);
        }
      // ------------------------------------


      this->_comm.barrier();


      // ------------------------------------
      // Solve for solution at integration points
      if( this->_comm.rank() == 0)
        std::cout << "     --> Solving at " << _nIntegrationPoints 
          << " integration points " << std::endl;
      
      std::vector<std::map<std::string,T_P> > myContribs ;
      std::vector<std::vector<double> > myPolyVals ;

      // solve for my integration points;
      for(unsigned int pt=0; pt < nPts; pt++)
      {
        if (AGNOS_DEBUG)
          std::cout << "DEBUG: beginning of pt" << std::endl;


        if (AGNOS_DEBUG)
          std::cout << "DEBUG: pt (pre computeContribution)" << std::endl;

        // compute the contribution
        myContribs.push_back(
          computeContribution( this->_physics, _integrationIndices[pt]) 
          );

        if (AGNOS_DEBUG)
          std::cout << "DEBUG: pt (post computeContribution)" << std::endl;

        myPolyVals.push_back( 
            evaluateBasis( this->_indexSet,
              _integrationPoints[_integrationIndices[pt]])
            );

        if (AGNOS_DEBUG)
          std::cout << "DEBUG: end of pt" << std::endl;
      } // end for pt 

      // clean up primary evals
      this->_primaryEvaluations.clear();

      this->_comm.barrier();


      // add my poly evaluations to global matrix
      if (AGNOS_DEBUG)
        std::cout << "DEBUG: begin copy poly evals" << std::endl;
      for(unsigned int i=0; i<myPolyVals.size(); i++)
        for(unsigned int j=0; j<myPolyVals[i].size(); j++)
        {
          polyValues->set(
              j, 
              _integrationIndices[i],
              myPolyVals[i][j] ) ;
        }
      polyValues->close();
      if (AGNOS_DEBUG)
        std::cout << "DEBUG: end copy poly evals" << std::endl;


      
      this->_comm.barrier();
      

      // need to know size of solution vector on all processes (in case some
      // don't have any work but will still be called in reduce operation)
      unsigned int tempSize;
      this->_solSize.clear();
      for (id=this->_solutionNames.begin();
          id!=this->_solutionNames.end(); id++)
      {
        if (this->_comm.rank() == 0)
          tempSize = myContribs[0][*id].size();

        if(AGNOS_DEBUG)
          std::cout << "DEBUG: surrogate build -- solSize[" << *id << "]:" <<
            tempSize << std::endl;

        this->_comm.broadcast(tempSize);

        if(AGNOS_DEBUG)
          std::cout << "DEBUG: surrogate build --(post_broadcast) solSize[" << *id << "]:" <<
            tempSize << std::endl;

        this->_solSize.insert( std::pair<std::string,unsigned int>(
              *id,tempSize ) );

      }



      // WAIT FOR ALL PROCESSES TO CATCH UP
      this->_comm.barrier();



      // --------------
      // compute my contribution to each coefficient
      // sum all processes contributions
      // save if its one of my coefficients
      // --------------
      
      if( this->_comm.rank() == 0)
        std::cout << "     --> Computing " << totalNCoeff << " coefficients" <<
          std::endl;


      //-- loop through sols
      for (id=this->_solutionNames.begin();
          id!=this->_solutionNames.end(); id++)
      {
        if(AGNOS_DEBUG)
          std::cout << "DEBUG: surrogate build computing coeffs: " << *id << std::endl ;

        unsigned int solSize = this->_solSize[*id] ;

        // get proper sizing to prevent PetscMatrix from being oversized after
        // multiplication
        // NOTE:  we will still fill all solComp for this procs integrationPts
        unsigned int minComp = solSize / this->_comm.size() ;
        unsigned int remComp = solSize % this->_comm.size() ;
        unsigned int nComp  =  minComp + ( this->_comm.rank() < remComp ) ;

        // Matrix to save solutions 
        std::shared_ptr<DistMatrix> solContrib( new DistMatrix(this->_comm) );
        solContrib->init(
            this->_nIntegrationPoints,  // global dim m
            solSize,                    // global dim n
            nPts,                       // local dim m
            nComp,                      // local dim n
            nComp,                      // # non-zero/row on proc
            solSize-nComp               // # non-zero/row off proc
            );
        solContrib->zero();

        // ensure that everthing is sized correctly
        if ( nPts > 0)
        {
          assert( solContrib->row_start() == this->_integrationIndices.front() );
          assert( solContrib->row_stop() == this->_integrationIndices.back()+1  );
        }

        // copy solution into solution matrix
        for (unsigned int i=0; i<myContribs.size(); i++)
          for (unsigned int j=0; j<solSize; j++)
            solContrib->set(
                _integrationIndices[i], 
                j,
                myContribs[i][*id](j) );
        solContrib->close();


        // Matrix product with plynomial values to compute coefficients
        Mat resultMat ;
        PetscErrorCode ierr;
        ierr = MatMatMult( 
            polyValues->mat(),
            solContrib->mat(),
            MAT_INITIAL_MATRIX,PETSC_DEFAULT,
            &resultMat
            );


        // Create libmesh DistMatrix from Petsc Mat 
        std::shared_ptr<DistMatrix> coeffMatrix(new
            DistMatrix(resultMat,this->_comm) ) ;

        /* std::cout << "coeffMatrx.m: " << coeffMatrix->m() << std::endl; */
        /* std::cout << "coeffMatrx.n: " << coeffMatrix->n() << std::endl; */
        /* for (unsigned int i=coeffMatrix->row_start(); i<coeffMatrix->row_stop(); */
        /*     i++) */
        /*   for (unsigned int j=0; j<coeffMatrix->n(); j++) */
        /*     std::cout << "coeff("<<i<<","<<j<<"): " << (*coeffMatrix)(i,j) */ 
        /*          << "rank: " << this->_comm.rank() << std::endl; */


        // Save solution coefficients in coefficient container
        this->_coefficients.insert( std::pair<std::string,std::shared_ptr<DistMatrix> >(
              *id, 
              coeffMatrix
              )
            );

      } // id



      if(AGNOS_DEBUG)
        std::cout << "DEBUG: leaving surrogate build " << std::endl ;

      return;
    } 

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, T_P >
    SurrogatePseudoSpectral<T_S,T_P>::computeContribution(
        std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
        unsigned int index
        )
    {

      if (AGNOS_DEBUG)
        std::cout << "DEBUG: computeContribution( ) beginning" << std::endl;
      

      // get data for this integration point
      std::map< std::string, T_P > contrib ;
      if ( ! this->_primaryEvaluations.empty() )
        contrib =
          this->_primaryEvaluations[index-this->_integrationIndices.front()];
      T_S integrationPoint = _integrationPoints[index];
      double integrationWeight = _integrationWeights[index];

      
      if (AGNOS_DEBUG)
        std::cout << "DEBUG: computeContribution( ) pre compute()" << std::endl;

      // get solution for this integration point
      physics->compute( integrationPoint, contrib );


      // clear out any evaluations that were provided by evalSurrogate, and that
      // aren't in solutionNames
      std::set<std::string>::iterator evalName = this->_evalNames.begin() ;
      for( ; evalName != this->_evalNames.end(); evalName++)
      {
        if ( !this->_solutionNames.count(*evalName) ) 
          contrib.erase( *evalName );
      }
      

      if (AGNOS_DEBUG)
      {
        std::cout << "DEBUG: computeContribution( ) post compute()" << std::endl;
      }

      // -- loop over solution names
      std::set<std::string>::iterator id = this->_solutionNames.begin();
      for(;id!=this->_solutionNames.end();++id)
        contrib[*id].scale( integrationWeight );


      if (AGNOS_DEBUG)
        std::cout << "DEBUG: computeContribution( ) ending" << std::endl;
      return contrib;
    }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S, class T_P>
    std::map<std::string, T_P> SurrogatePseudoSpectral<T_S,T_P>::evaluate( 
        std::set< std::string > solutionNames,
        T_S& parameterValues, /**< parameter values to evaluate*/
        bool saveLocal  /**< save solution locally after evaluation*/
        ) const
    {
      if(AGNOS_DEBUG)
        std::cout << "DEBUG: entering surrogate evaluate" << std::endl;

      // initalize some data structure
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::map< std::string, T_P> surrogateValue;

      // initalize polyValues data structure
      // NOTE: we will fill all polys on each proc, 
      //        i.e. global m = comm.size
      //        local m = 1
      DistMatrix polyValues(this->_comm);
      polyValues.init( 
          this->_comm.size(),           // global dim m
          totalNCoeff,                  // global dim n
          saveLocal,                    // local dim m
          this->_coeffIndices.size(),   // local dim n
          totalNCoeff,                  // # non-zero/row on proc
          totalNCoeff                   // # non-zero/row off proc
          );
      polyValues.zero();



      if (saveLocal)
      {
        // make sure we only have one local row
        assert( (polyValues.row_stop() - polyValues.row_start()) == 1);
        
        // evaluate all basis polys at given parameterValue
        std::vector<double> myPolyVals = evaluateBasis( 
            this->_indexSet,
            parameterValues 
            );
      
        // copy poly vals into Petsc data structure
        for(unsigned int i=0; i<myPolyVals.size(); i++)
          polyValues.set( polyValues.row_start(), i , myPolyVals[i] ) ;
      }
      else
      {
        // make sure we have NO local rows
        assert( (polyValues.row_stop() - polyValues.row_start()) == 0);
      }
      polyValues.close();



      this->_comm.barrier();



      // loop over all solutions requested
      std::set<std::string>::iterator id = solutionNames.begin();
      for (; id != solutionNames.end(); id++)
      {
        // reference for solution size
        unsigned int solSize 
          = const_cast<SurrogatePseudoSpectral<T_S,T_P>*>(this)->_solSize[*id] ;

        if(AGNOS_DEBUG)
          std::cout << "DEBUG: evaluating surrogate model for: " << *id <<
            "with size:" <<  solSize << std::endl;


        // compute matrix matrix product b/w polyValues and coefficient matrix
        Mat resultMat ;
        PetscErrorCode ierr;
        ierr = MatMatMult( 
            polyValues.mat(),
            (const_cast<SurrogatePseudoSpectral<T_S,T_P>*>(this)->_coefficients[*id])->mat(),
            MAT_INITIAL_MATRIX,PETSC_DEFAULT,
            &resultMat
            );

        // initiate matrix from MM product result
        std::shared_ptr<DistMatrix> coeffMatrix(
            new DistMatrix(resultMat,this->_comm) ) ;

        // localize result on each processor
        if (saveLocal)
        {
          std::vector<double> localResult;
          for (unsigned int i=coeffMatrix->row_start(); 
              i<coeffMatrix->row_stop(); i++)
          {
            localResult.resize(solSize);
            for (unsigned int j=0; j<solSize; j++)
              localResult[j]=(*coeffMatrix)(i,j)  ;
          }
          
          // save result in return map
          surrogateValue.insert(
              std::pair< std::string, T_P >(*id, T_P(localResult) ) 
              ) ;

        }


      } // solNames

      if(AGNOS_DEBUG)
        std::cout << "DEBUG: leaving surrogate evaluate" << std::endl;
      return surrogateValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationPoints( 
      std::ostream& out ) const
  {
    if (this->_comm.rank() == 0)
    {
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Integration points " << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < _nIntegrationPoints ; ix++)  
      {
        for(int iy=0; iy < this->_dimension; iy++)
        {
          out << std::scientific << std::setprecision(5) << std::setw(12)
            << _integrationPoints[ix](iy) << "  ";
        }
        out << std::endl;
      }
      out << "#----------------------------------------------------" <<
        std::endl;
    }

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::prettyPrintIntegrationPoints( ) const
  {
    if (this->_comm.rank() == 0)
    {
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Integration points " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      std::cout << "   \\ x        " << std::endl;
      std::cout << "    \\" ;
      for(unsigned int dim=0; dim < this->_dimension; dim++)
        std::cout << std::setw(12) << "x_" << dim << " " ;
      std::cout << std::endl;
      std::cout << "  id \\  " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < _nIntegrationPoints ; ix++)  
      {
        std::cout << std::setw(5) << ix << " |   ";
        for(int iy=0; iy < this->_dimension; iy++)
        {
          std::cout << std::scientific << std::setprecision(5) << std::setw(12)
            << _integrationPoints[ix](iy) << "  ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    return;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationWeights( 
      std::ostream& out ) const
  {
    if (this->_comm.rank() == 0)
    {
      double sum = 0.0;
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Integration weights " << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
      for(int ip=0; ip < _nIntegrationPoints; ip++){  
        out << _integrationWeights[ip] << "  ";
        out << std::endl;
        sum += _integrationWeights[ip];
      }

      out << std::endl;
      out << "#Sum = " << sum << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
    }
    return;
  }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
    template<class T_S,class T_P>
    void SurrogatePseudoSpectral<T_S,T_P>::prettyPrintIntegrationWeights( ) const
    {
      if (this->_comm.rank() == 0)
      {
        double sum = 0.0;
        std::cout << std::endl;
        std::cout << "====================================================" <<
          std::endl;
        std::cout << " Integration weights " << std::endl;
        std::cout << "----------------------------------------------------" <<
          std::endl;
        for(int ip=0; ip < _nIntegrationPoints; ip++){  
          std::cout << std::setw(5) << ip << "   |   ";
          std::cout << _integrationWeights[ip] << "  ";
          std::cout << std::endl;
          sum += _integrationWeights[ip];
        }

        std::cout << std::endl;
        std::cout << "Sum = " << sum << std::endl;
        std::cout << std::endl;
      }
      return;
    }


  /********************************************//**
   * \brief 
   ***********************************************/
    template<class T_S, class T_P> 
      void SurrogatePseudoSpectral<T_S,T_P>::printIndexSet( std::ostream& out ) const
      {
        if (this->_comm.rank() == 0)
        {
          out << std::endl;
          out << "#====================================================" <<
            std::endl;
          out << "# Index Set" << std::endl;
          out << "#----------------------------------------------------" <<
          std::endl;

          for (unsigned int i=0; i< _indexSet.size(); i++)
          {
            for (unsigned int j=0; j< _indexSet[i].size(); j++)
            {
              out << std::setw(5) << _indexSet[i][j] << " " ;
            }
            out << std::endl;
          }
          out << "#----------------------------------------------------" <<
            std::endl;
        }

      return ;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P> 
    void SurrogatePseudoSpectral<T_S,T_P>::prettyPrintIndexSet( ) const
    {
      if (this->_comm.rank() == 0)
      {
        std::cout << std::endl;
        std::cout << "====================================================" <<
          std::endl;
        std::cout << " Index Set" << std::endl;
        std::cout << "----------------------------------------------------" <<
        std::endl;
      std::cout << "   \\ dir      " << std::endl;
      std::cout << "    \\        " ;
      for(unsigned int dim=0; dim < this->_dimension; dim++)
        std::cout << std::setw(4) << "xi_" << dim << " " ;
      std::cout << std::endl;
      std::cout << "  id \\  " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
        for (unsigned int i=0; i< _indexSet.size(); i++)
        {
        std::cout << std::setw(5) << i << "   |   ";
          for (unsigned int j=0; j< _indexSet[i].size(); j++)
          {
            std::cout << std::setw(5) << _indexSet[i][j] << " " ;
          }
          std::cout << std::endl;
        }

      }
      return ;
    }

  /********************************************//**
   * \brief calculate l2 norm by integration
   ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, T_P> SurrogatePseudoSpectral<T_S,T_P>::l2Norm(
        std::set< std::string > solutionNames
        )
    {

      // initialize to zero
      std::map< std::string, T_P> l2norm;
      std::set<std::string>::iterator id = solutionNames.begin();
      for (; id != solutionNames.end(); id++)
      {
        l2norm.insert( std::pair<std::string,T_P>(
              *id, T_P( std::vector<double>(this->_solSize[*id],0.) ) )
            );
      }

      //---------
      //get integration points for higher order quad rule
      std::vector<unsigned int> integrationOrder = this->_order;
      for(unsigned int i=0; i<integrationOrder.size();i++)
        integrationOrder[i] += 2;

      QuadratureTensorProduct integrationQuadratureRule(
          this->_parameters, integrationOrder );

      unsigned int dimension = this->_dimension;
      unsigned int nQuadPoints = integrationQuadratureRule.getNQuadPoints();
      double** quadPoints = integrationQuadratureRule.getQuadPoints();
      double* quadWeights = integrationQuadratureRule.getQuadWeights();



      //----------
      // loop over quad points
      for(unsigned int i=0; i<nQuadPoints; i++)
      {
        T_S paramValues(dimension);
        for(unsigned int j=0; j<paramValues.size(); j++)
          paramValues(j) = quadPoints[i][j];

        std::map<std::string,T_P> pointValue 
          = this->evaluate(solutionNames,paramValues);
        for (id=solutionNames.begin(); id != solutionNames.end(); id++)
        {
          for(unsigned int j=0; j<l2norm[*id].size();j++)
            l2norm[*id](j) += pointValue[*id](j) * pointValue[*id](j) 
              * quadWeights[i] ; 
        }

      } // end loop over quad points

      
      // take square root
      for (id=solutionNames.begin(); id != solutionNames.end(); id++)
      {

        for(unsigned int j=0; j<l2norm[*id].size();j++)
          l2norm[*id](j) =  std::sqrt( l2norm[*id](j) );
      }


      return l2norm;
    }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    double SurrogatePseudoSpectral<T_S,T_P>::l2NormDifference(
        SurrogateModel<T_S,T_P>& comparisonModel,
        std::string solutionName )
    {

      // make sure solutionName is present in both models
      // evaluate at integration points
      //
      //
      // initialize to zero
      double l2norm = 0;

      if (this->_solutionNames.count(solutionName) == 0) 
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " base model of l2NormDifference( ... ). \n "
          << std::endl;
        exit(1);
      }
      if (comparisonModel.getSolutionNames().count(solutionName) == 0 )
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " comparison model of l2NormDifference( ... ). \n "
          << std::endl;
        exit(1);
      }

      //---------
      //get integration points for higher order quad rule
      std::vector<unsigned int> integrationOrder = this->_order;
      std::vector<unsigned int> comparisonOrder 
        = comparisonModel.getExpansionOrder() ;
      for(unsigned int i=0; i<integrationOrder.size();i++)
        integrationOrder[i] 
          = std::max(integrationOrder[i],comparisonOrder[i]) 
          + 2;

      QuadratureTensorProduct integrationQuadratureRule(
          this->_parameters, integrationOrder );

      unsigned int dimension = this->_dimension;
      unsigned int nQuadPoints = integrationQuadratureRule.getNQuadPoints();
      double** quadPoints = integrationQuadratureRule.getQuadPoints();
      double* quadWeights = integrationQuadratureRule.getQuadWeights();



      //----------
      // loop over quad points
      for(unsigned int i=0; i<nQuadPoints; i++)
      {
        T_S paramValues(dimension);
        for(unsigned int j=0; j<paramValues.size(); j++)
          paramValues(j) = quadPoints[i][j];

        T_P diffVec = this->evaluate(solutionName,paramValues);
        diffVec -= comparisonModel.evaluate(solutionName,paramValues);


        l2norm += diffVec.dot( diffVec )  * quadWeights[i] ; 
        /* l2norm += 1  * quadWeights[i] ; */ 

      } // end loop over quad points

      
      // take square root
      l2norm =  std::sqrt( l2norm );


      return l2norm;

    }
  
}
#endif // SURROGATE_PSEUDO_SPECTRAL_H


