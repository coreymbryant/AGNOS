
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
        std::shared_ptr<SurrogateModel<T_S,T_P> > primarySurrogate, 
        std::vector<unsigned int> increaseOrder ,
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
    SurrogatePseudoSpectral<T_S,T_P>::indexSet( ) const
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
      // First call initialize to set up data structures
      this->initialize();
      // This is separated from the routine that actually computes contribution
      // so that we can group surrogate models together later and wll that needs
      // to be defined is build routine based on computeContribution( )
      //
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::set< std::string >::iterator id;
      
      // ------------------------------------
      // Determine processes position in the MPI groups
      int globalRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
      // ------------------------------------
      

      // output progress
      if( this->_groupRank == 0)
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
      unsigned int minPts = _nIntegrationPoints / this->_nPhysicsGroups ;
      unsigned int remPts = _nIntegrationPoints % this->_nPhysicsGroups ;
      unsigned int nPts  =  minPts + ( this->_physicsGroup < remPts ) ;
      unsigned int intPtsStart = this->_physicsGroup*(minPts) 
        + std::min( (unsigned int)this->_physicsGroup, remPts ) ;

      this->_integrationIndices.clear();
      for (unsigned int i=0; i < nPts; i++)
        _integrationIndices.push_back( intPtsStart + i );

      assert( nPts == _integrationIndices.size() );

      // assign the appropriate coeff to each processor
      unsigned int minCoeff = totalNCoeff / this->_nPhysicsGroups ;
      unsigned int remCoeff = totalNCoeff % this->_nPhysicsGroups ;
      unsigned int nCoeffs  =  minCoeff + ( this->_physicsGroup < remCoeff ) ;
      unsigned int coeffStart = this->_physicsGroup*(minCoeff) 
        + std::min( (unsigned int)this->_physicsGroup, remCoeff ) ;

      this->_coeffIndices.clear();
      for (unsigned int i=0; i < nCoeffs; i++)
        this->_coeffIndices.push_back( coeffStart + i );
      // ------------------------------------


      /* if (AGNOS_DEBUG) */
        std::cout << "physicsGroup: " << this->_physicsGroup
          << " groupRank: " << this->_groupRank 
          << " global rank: " << globalRank 
          << " nPts: " << this->_integrationIndices.size()
          << " comm.size: " << this->_nPhysicsGroups 
          << std::endl;



      // ------------------------------------
      // Initialize data structures
      // TODO: issues start here if libMesh not built with
      // --disable-default-commWorld
      std::shared_ptr<DistMatrix> polyValues;
      if (this->_groupRank ==0)
      {
        polyValues.reset(new DistMatrix(this->_comm));
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
      }

      // ------------------------------------

      
      /* this->_comm.barrier(); */


      // ------------------------------------
      // Solve for solution at integration points
      if( this->_groupRank == 0)
        std::cout << "     --> Solving at " << _nIntegrationPoints 
          << " integration points " << std::endl;
      
      // ------------------------------------
      // Get primalSurrogate evaluations if needed
      // ------------------------------------
      this->_primaryEvaluations.clear();
      this->_primaryEvaluations.reserve( nPts );
      if ( ! this->_evalNames.empty() )
      {
        // loop through all integration points
        for(unsigned int pt=0; pt < this->_nIntegrationPoints; pt++)
        {
          if (AGNOS_DEBUG)
            std::cout << "DEBUG: evaluating primarySurrogate at pt: " << pt 
              << std::endl;
          
          // determine if its one of this processes points
          bool localPoint = false;
          if (nPts > 0)
            localPoint = ( 
              (pt >= this->_integrationIndices.front()) 
              && 
              (pt <= this->_integrationIndices.back()) 
              ) ;

          std::map<std::string,T_P> evaluations;
          if (this->_groupRank==0)
            evaluations  = this->_evalSurrogate->evaluate(
                this->_evalNames, _integrationPoints[pt], localPoint  ) ;
          if (localPoint)
            this->_primaryEvaluations.push_back(evaluations);
        }
      }
      // ------------------------------------

      // ------------------------------------
      // TODO update since its not necessary now
      // loop through all integration points, even if not mine so that
      // evaluations of primary surrogate operate collectively.
      std::map<std::string, std::shared_ptr<DistMatrix> > solContrib;
      std::map<std::string,T_P> myContribs ;
      for(unsigned int pt=0; pt < this->_nIntegrationPoints; pt++)
      {
        if ( pt%10 == 0)
          std::cout << "pt " << pt << std::endl;

        if (AGNOS_DEBUG)
          std::cout << "DEBUG: beginning of pt" << std::endl;

        bool localPoint = false;
        if (nPts > 0)
          localPoint = ( 
            (pt >= this->_integrationIndices.front()) 
            && 
            (pt <= this->_integrationIndices.back()) 
            ) ;
        
      
       
        
        /* this->_comm.barrier(); */


        // ------------------------------------
        // compute the first coefficient and broadcast size to all processes
        // (in case some don't have any work but will still be called in reduce
        // operation)
        //  - this is also when we will initialize the solContrib DistMatrix
        if (pt == 0)
        {
          if (this->_physicsGroup == 0)
          {
            myContribs = computeContribution( this->_physics, 0) ;
          }


          // let proc 0 catch up
          /* this->_comm.barrier(); */

          // Matrix to save solutions 
          unsigned int tempSize;
          this->_solSize.clear();
          for (id=this->_solutionNames.begin();
              id!=this->_solutionNames.end(); id++)
          {
            if (this->_physicsGroup == 0)
              tempSize = myContribs[*id].size();

            if(AGNOS_DEBUG)
              std::cout << "DEBUG: surrogate build -- solSize[" << *id << "]:" <<
                tempSize << std::endl;

            this->_comm.broadcast(tempSize);

            if(AGNOS_DEBUG)
              std::cout << "DEBUG: surrogate build --(post_broadcast) solSize["
                << *id << "]:" << tempSize << std::endl;

            this->_solSize.insert( std::pair<std::string,unsigned int>(
                  *id,tempSize ) );

            
            // add a matrix for this solution and initalize to correct size
            if (this->_groupRank == 0)
            {
              solContrib.insert( 
                  std::pair<std::string, std::shared_ptr<DistMatrix> >(
                    *id, std::shared_ptr<DistMatrix>(new DistMatrix(this->_comm))  
                    ) 
                  );
            }

            unsigned int solSize = this->_solSize[*id] ;

            // get proper sizing to prevent PetscMatrix from being oversized after
            // multiplication
            // NOTE:  we will still fill all solComp for this procs integrationPts
            unsigned int minComp = solSize / this->_nPhysicsGroups ;
            unsigned int remComp = solSize % this->_nPhysicsGroups ;
            unsigned int nComp  =  minComp + ( this->_physicsGroup < remComp ) ;

            if (this->_groupRank == 0)
            {
              // Matrix to save solutions 
              solContrib[*id]->init(
                  this->_nIntegrationPoints,  // global dim m
                  solSize,                    // global dim n
                  nPts,                       // local dim m
                  nComp,                      // local dim n
                  nComp,                      // # non-zero/row on proc
                  solSize-nComp               // # non-zero/row off proc
                  );
              solContrib[*id]->zero();

              // ensure that everthing is sized correctly
              if ( nPts > 0)
              {
                assert( solContrib[*id]->row_start() == this->_integrationIndices.front() );
                assert( solContrib[*id]->row_stop() == this->_integrationIndices.back()+1  );
              }

              // copy first coefficient into matrix
              if (this->_physicsGroup == 0)
              {
                for (unsigned int j=0; j<solSize; j++)
                  solContrib[*id]->set(
                      0, 
                      j,
                      myContribs[*id](j) );
              }
            }
          } // end of solution names
        
        } // end of if pt==0

        // otherwise solve if its one of my points
        else if ( localPoint )
        {
          if (AGNOS_DEBUG)
            std::cout << "DEBUG: pt (pre computeContribution)" << std::endl;

          // compute the contribution
          std::map<std::string,T_P> myContribs =
            computeContribution( this->_physics, pt) ;

          if (this->_groupRank==0)
          {
            // and save into solution matrix
            for (id=this->_solutionNames.begin();
                id!=this->_solutionNames.end(); id++)
            {
              for (unsigned int j=0; j<this->_solSize[*id]; j++)
                solContrib[*id]->set(
                    pt, 
                    j,
                    myContribs[*id](j) );
            }
          }

          if (AGNOS_DEBUG)
            std::cout << "DEBUG: pt (post computeContribution)" << std::endl;
        } // end of if local point and not pt==0



        // now compute poly values
        if (localPoint)
        {
          if (AGNOS_DEBUG)
            std::cout << "DEBUG: begin copy poly evals" << std::endl;
          if(this->_groupRank==0)
          {
            std::vector<double> myPolyVals =  
              evaluateBasis( this->_indexSet,
                  _integrationPoints[pt]) ;
            for(unsigned int j=0; j<myPolyVals.size(); j++)
              polyValues->set(
                  j, 
                  pt,
                  myPolyVals[j] ) ;
          }
          if (AGNOS_DEBUG)
            std::cout << "DEBUG: end copy poly evals" << std::endl;
        }


        /* this->_comm.barrier(); */


        if (AGNOS_DEBUG)
          std::cout << "DEBUG: end of pt" << std::endl;
      }// end of pts


      // clean up primary evals
      this->_primaryEvaluations.clear();


      /* this->_comm.barrier(); */
      if (this->_groupRank==0)
        polyValues->close();


      // --------------
      // compute my contribution to each coefficient
      // sum all processes contributions
      // save if its one of my coefficients
      // --------------
      
      if( this->_groupRank == 0)
        std::cout << "     --> Computing " << totalNCoeff << " coefficients" <<
          std::endl;


      if (this->_groupRank==0)
      {
        //-- loop through sols
        for (id=this->_solutionNames.begin();
            id!=this->_solutionNames.end(); id++)
        {
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: surrogate build computing coeffs: " << *id << std::endl ;

          // close the matrix now, we are done adding to it
          solContrib[*id]->close();

          // Matrix product with plynomial values to compute coefficients
          Mat resultMat ;
          PetscErrorCode ierr;
          ierr = MatMatMult( 
              polyValues->mat(),
              solContrib[*id]->mat(),
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
      }



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

      PetscErrorCode ierr;

      if(this->_groupRank==0)
      {

        // initalize polyValues data structure
        Vec polyValuesVec;
        ierr = VecCreateMPI(
            this->_comm.get(),  // comm
            PETSC_DECIDE,       // local size (let petsc decide)
            totalNCoeff,        // global size
            &polyValuesVec);    // vector
        Vector polyValues(polyValuesVec,this->_comm);
        polyValues.zero();

          
        // evaluate all basis polys at given parameterValue
        std::vector<double> myPolyVals = evaluateBasis( 
            this->_indexSet,
            parameterValues 
            );
        
        // copy poly vals into Petsc data structure
        for(unsigned int i=0; i<polyValues.local_size(); i++)
          polyValues.set( i+polyValues.first_local_index() ,
              myPolyVals[i+polyValues.first_local_index()] ) ;

        // close vector
        polyValues.close();


        // loop over all solutions requested
        std::set<std::string>::iterator id = solutionNames.begin();
        for (; id != solutionNames.end(); id++)
        {
          // Make sure we have all the requested solution names
          if
            (!const_cast<SurrogatePseudoSpectral<T_S,T_P>*>(this)->getSolutionNames().count(*id))
          {
            std::cout << std::endl;
            std::cerr << 
              " ERROR: requested evaluation for solution that isn't present  "
              << std::endl;
            std::cout << std::endl;
            std::abort();
          }

          // reference for solution size
          unsigned int solSize 
            = const_cast<SurrogatePseudoSpectral<T_S,T_P>*>(this)->_solSize[*id] ;


          if(AGNOS_DEBUG)
            std::cout << "DEBUG: evaluating surrogate model for: " << *id <<
              "with size:" <<  solSize << std::endl;


          // create vector to store result of evaluation
          Vec resultVec;;
          ierr = VecCreateMPI(
              this->_comm.get(), 
              PETSC_DECIDE,
              solSize,
              &resultVec);
          Vector result(resultVec,this->_comm);
          result.zero();

          // compute  matrix vector product b/w polyValues and coefficient matrix
          ierr = MatMultTranspose( 
              (const_cast<SurrogatePseudoSpectral<T_S,T_P>*>(this)->_coefficients[*id])->mat(),
              polyValues.vec(),
              result.vec()
              );


          // localize result on each processor
          std::vector<double> localResult(solSize,0.);
          result.localize(localResult);
            
          if (saveLocal)
          {
            // save result in return map
            surrogateValue.insert(
                std::pair< std::string, T_P >(*id, T_P(localResult) ) 
                ) ;

          }

        } // solNames
      } // if groupRank ==0

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
    if (this->_groupRank == 0)
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
    if (this->_groupRank == 0)
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
    if (this->_groupRank == 0)
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
      if (this->_groupRank == 0)
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
        if (this->_groupRank == 0)
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
      if (this->_groupRank == 0)
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

      if(this->_groupRank==0)
      {
        //---------
        //get integration points for higher order quad rule
        std::vector<unsigned int> integrationOrder = this->_order;
        for(unsigned int i=0; i<integrationOrder.size();i++)
          integrationOrder[i] += 0;

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

      } // end groupRank == 0

      
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

      if(this->_groupRank==0)
      {
        //---------
        //get integration points for higher order quad rule
        std::vector<unsigned int> integrationOrder = this->_order;
        std::vector<unsigned int> comparisonOrder 
          = comparisonModel.getExpansionOrder() ;
        for(unsigned int i=0; i<integrationOrder.size();i++)
          integrationOrder[i] 
            = std::max(integrationOrder[i],comparisonOrder[i]) 
            + 0;

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
          /* std::cout << " test: base model eval " << std::endl; */
          diffVec -= comparisonModel.evaluate(solutionName,paramValues);
          /* std::cout << " test: comparison model eval " << std::endl; */


          l2norm += diffVec.dot( diffVec )  * quadWeights[i] ; 
          /* l2norm += 1  * quadWeights[i] ; */ 

        } // end loop over quad points

      } // end of if groupRank == 0

      
      // take square root
      l2norm =  std::sqrt( l2norm );


      return l2norm;

    }
  
}
#endif // SURROGATE_PSEUDO_SPECTRAL_H


