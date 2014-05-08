
#include "EvaluatorPseudoSpectral.h"

namespace AGNOS
{
  /********************************************//**
   * \brief Constructor
   ***********************************************/
  template<class T_S,class T_P>
    EvaluatorPseudoSpectral<T_S,T_P>::EvaluatorPseudoSpectral(
        const Communicator&               comm,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&          order,
        std::set<std::string> computeSolutions 
        )
    : SurrogateModelBase<T_S,T_P>(comm,parameters,order,computeSolutions),
      SurrogateEvaluator<T_S,T_P>(comm,parameters,order,computeSolutions)
    {
    }

  /********************************************//**
   * \brief Destructor
   ***********************************************/
  template<class T_S,class T_P>
    EvaluatorPseudoSpectral<T_S,T_P>::~EvaluatorPseudoSpectral()
    { 
    }

  /********************************************//**
   * \brief Build routine. Sets indexSet and Coefficients 
   ***********************************************/
  template<class T_S,class T_P>
    void EvaluatorPseudoSpectral<T_S,T_P>::build( 
        std::vector< std::vector<unsigned int> > indexSet,
        std::map< std::string, LocalMatrix >  coefficients
        ) 
    {

      // set up some of the members that aren't going to be set explicitly here
      this->_totalNCoeff = coefficients.begin()->second.m(); 

      // set up a map of distCoefficients to feed to other build routine
      std::map<std::string, std::shared_ptr<DistMatrix> > distCoefficients;

      // get my rank
      unsigned int myRank = this->_physicsGroup;

      // if rank 0 in my physics group
      if (this->_groupRank==0)
      {

        // loop over all solution names
        std::set<std::string>::iterator id = this->_solutionNames.begin();
        for (; id!=this->_solutionNames.end(); ++id)
        {
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: build routine: id: " << *id << std::endl;

          distCoefficients.insert(
              std::pair<std::string,std::shared_ptr<DistMatrix> >(
                *id, std::shared_ptr<DistMatrix>(new DistMatrix(this->_comm))  
                )
              ) ;

          this->_solSize[*id] = coefficients[*id].n();

          unsigned int solSize = this->_solSize[*id];
          unsigned int totalNCoeff = this->_totalNCoeff;
          unsigned int minCoeff = totalNCoeff / this->_nPhysicsGroups ;
          unsigned int remCoeff = totalNCoeff % this->_nPhysicsGroups ;
          unsigned int nCoeffs  =  minCoeff + ( this->_physicsGroup < remCoeff ) ;
          unsigned int coeffStart = this->_physicsGroup*(minCoeff) 
            + std::min( (unsigned int)this->_physicsGroup, remCoeff ) ;
          unsigned int minComp = solSize / this->_nPhysicsGroups ;
          unsigned int remComp = solSize % this->_nPhysicsGroups ;
          unsigned int nComp  =  minComp + ( this->_physicsGroup < remComp ) ;

          this->_coeffIndices.clear();
          for (unsigned int i=0; i < nCoeffs; i++)
            this->_coeffIndices.push_back( coeffStart + i );


          distCoefficients[*id]->init(
              this->_totalNCoeff, // global dim M
              solSize,            // global dim N
              nCoeffs,            // local dim m
              nComp,              // local dim n
              nComp,              // # non-zero/row on proc
              solSize-nComp       // # non-zero/row off proc
              );
          distCoefficients[*id]->zero();
          // ensure that everthing is sized correctly
          if ( nCoeffs > 0)
          {
            assert( distCoefficients[*id]->row_start() == this->_coeffIndices.front() );
            assert( distCoefficients[*id]->row_stop() == this->_coeffIndices.back()+1  );
          }

          if (myRank == 0)
          {
            /* std::cout << "post myrank = 0 " << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.front(): " << */
            /*   this->_coeffIndices.front() << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.back(): " << */
            /*   this->_coeffIndices.back() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.row_start(): " << */
            /*   distCoefficients[*id]->row_start() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.row_stop(): " << */
            /*   distCoefficients[*id]->row_stop() << std::endl; */
            
            // save my coeffs
            for (unsigned int i=this->_coeffIndices.front(); 
                  i<this->_coeffIndices.back()+1; i++)
              for (unsigned int j=0; j<solSize; j++)
              {
                /* std::cout<<"solCoeff("<<i<<","<<j<<"):" */
                /*   << (*this->_coefficients[*id])(i,j) << std::endl; */
                distCoefficients[*id]->set(i,j, coefficients[*id](i,j))  ;
              }
            distCoefficients[*id]->close();


            // Send coefficients to other procs
            for (unsigned int p=1; p<this->_comm.size(); p++)
            {
              std::cout << " rank 0 sending to rank " << p << std::endl;
              std::vector<unsigned int> jbuf;
              std::vector<double> cbuf;

              // first get coeffIndices for this proc
              this->_comm.receive(p, jbuf);

              /* std::cout << "rank 0: recv from rank: " << p << " jbuf.size: " << */
              /*   jbuf.size()<< std::endl; */


              // if proc p has non-empty index set send its coefficients
              if (!jbuf.empty())
              {
                for (unsigned int i=jbuf.front(); i<jbuf.back()+1; i++)
                  for (unsigned int j=0; j<solSize; j++)
                  {
                    /* std::cout << "rank " << myRank << " _coefficients(i,j): " << */
                    /*   (*this->_coefficients[*id])(i,j) << std::endl; */
                    cbuf.push_back( (coefficients[*id])(i,j) )  ;
                  }
              } // if jbuf not empty

              /* std::cout << "rank 0 cbuf.size: " << cbuf.size() << std::endl; */
              this->_comm.send(p,cbuf);

            } // for each proc p
          } // if rank 0
          else // receive my coeffs from proc 0
          {
            /* std::cout << "rank " << myRank << std::endl; */
            std::vector<unsigned int> jbuf;
            std::vector<double> cbuf;


            /* std::cout << " rank: " << myRank << " _coeffIndices.size: " << */
            /*   this->_coeffIndices.size() << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.front(): " << */
            /*   this->_coeffIndices.front() << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.back(): " << */
            /*   this->_coeffIndices.back() << std::endl; */

            // send my indices to proc 0
            if ( !this->_coeffIndices.empty() )
            {
              for (unsigned int i=this->_coeffIndices.front(); 
                    i<this->_coeffIndices.back()+1; i++)
              {
                jbuf.push_back(i) ;
              }
            } // if coeffIndices not empty

            /* std::cout << "rank " << myRank */ 
            /*   << " jbuf.size: " << jbuf.size() << std::endl; */
            this->_comm.send(0,jbuf);


            /* std::cout << "rank " << myRank << " recv from rank 0" << std::endl; */

            // receive my coeffs
            this->_comm.receive(0, cbuf);

            /* std::cout << "rank " << myRank << " recv from rank 0  cbuf.size: " */ 
            /*   << cbuf.size()<< std::endl; */


            if (!jbuf.empty())
            {
              for (unsigned int i=jbuf.front(); i<jbuf.back()+1; i++)
                for (unsigned int j=0; j<solSize; j++)
                {
                  unsigned int index = 
                    (i-jbuf.front()) * solSize  // due to j index
                    + j ;
                  /* std::cout<<"solCoeff("<<i<<","<<j<<"):" */
                  /*   << cbuf[index] << std::endl; */
                  distCoefficients[*id]->set(i,j, cbuf[index] ) ;
                } // loop over sol index
            } // if jbuf not empty

            distCoefficients[*id]->close();

          } // if not proc 0
        } // end id
      }

      this->build( indexSet, distCoefficients ) ;
      return;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P> 
    std::vector<double> EvaluatorPseudoSpectral<T_S,T_P>::evaluateBasis( 
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
 * \brief Initialization routine
 ***********************************************/
  template<class T_S, class T_P>
    std::map<std::string, T_P> EvaluatorPseudoSpectral<T_S,T_P>::evaluate( 
        std::set< std::string >  solutionNames,
        T_S&                        parameterValues,
        bool saveLocal ) const 
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
            (!this->getSolutionNames().count(*id))
          {
            std::cout << std::endl;
            std::cerr 
              << " ERROR: requested evaluation for solution " 
              << *id << ", " 
              << "which isn't present."
              << std::endl;
            std::cout << std::endl;
            std::abort();
          }

          // reference for solution size
          unsigned int solSize = 
            const_cast<EvaluatorPseudoSpectral<T_S,T_P>*>(this)->_solSize[*id] ;


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
              (const_cast<EvaluatorPseudoSpectral<T_S,T_P>*>(this)->_coefficients[*id])->mat(),
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

  template class 
    EvaluatorPseudoSpectral<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
} // namespace AGNOS
