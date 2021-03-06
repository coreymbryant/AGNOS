
#include "SurrogatePseudoSpectral.h"


namespace AGNOS
{
/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>& order,
        std::set<std::string> computeSolutions
        )
      : 
        SurrogateModelBase<T_S,T_P>(comm,parameters,order,computeSolutions),
        EvaluatorPseudoSpectral<T_S,T_P>(comm,parameters,order,computeSolutions),
        SurrogateModel<T_S,T_P>(comm,physics,parameters,order,computeSolutions)
    {
    }


/********************************************//*
 * \brief Secondary constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        std::shared_ptr<SurrogateModelBase<T_S,T_P> > primarySurrogate, 
        std::vector<unsigned int> increaseOrder ,
        unsigned int multiplyOrder ,
        std::set<std::string> evaluateSolutions,
        std::set<std::string> computeSolutions
        )
      : 
        SurrogateModelBase<T_S,T_P>(
            primarySurrogate->getComm(),
            primarySurrogate->getParameters(),
            primarySurrogate->getExpansionOrder(),
            computeSolutions),
        EvaluatorPseudoSpectral<T_S,T_P>(
            primarySurrogate->getComm(),
            primarySurrogate->getParameters(),
            primarySurrogate->getExpansionOrder(), 
            computeSolutions),
        SurrogateModel<T_S,T_P>(primarySurrogate, increaseOrder, multiplyOrder,
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
              this->evaluateBasis( this->_indexSet,
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

          /* std::cout << "polyValues.m: " << polyValues->m() << std::endl; */
          /* std::cout << "polyValues.n: " << polyValues->n() << std::endl; */
          /* for (unsigned int i=polyValues->row_start(); i<polyValues->row_stop(); */
          /*     i++) */
          /*   for (unsigned int j=0; j<polyValues->n(); j++) */
          /*     std::cout << "poly("<<i<<","<<j<<"): " << (*polyValues)(i,j) */ 
          /*          << "rank: " << this->_comm.rank() << std::endl; */


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
          /* std::cout << "solContrib[*id].m: " << solContrib[*id]->m() << std::endl; */
          /* std::cout << "solContrib[*id].n: " << solContrib[*id]->n() << std::endl; */
          /* for (unsigned int i=solContrib[*id]->row_start(); i<solContrib[*id]->row_stop(); */
          /*     i++) */
          /*   for (unsigned int j=0; j<solContrib[*id]->n(); j++) */
          /*     std::cout << "sol("<<i<<","<<j<<"): " << (*solContrib[*id])(i,j) */ 
          /*          << "rank: " << this->_comm.rank() << std::endl; */

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
      physics->compute( this->_solutionNames, integrationPoint, contrib );


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
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationPoints( 
      std::ostream& out ) const
  {
    if (this->_groupRank == 0)
    {
      out << std::endl;
      out << "#----------------------------------------------------" <<
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
      out << "#----------------------------------------------------" <<
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


  template class 
    SurrogatePseudoSpectral<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
}


