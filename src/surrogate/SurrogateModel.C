
#include "SurrogateModel.h"

namespace AGNOS
{


/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >           physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&          order,
        std::set<std::string> computeSolutions 
        )
      : 
        _physics(physics),
        SurrogateModelBase<T_S,T_P>(comm,parameters,order,computeSolutions)
    {
      int globalRank,globalSize;
      MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
      MPI_Comm_size(MPI_COMM_WORLD,&globalSize);
      this->_physicsGroup = globalRank / ( _physics->comm().size() ) ;
      this->_nPhysicsGroups = globalSize / ( _physics->comm().size() ) ;
      MPI_Comm_rank(_physics->comm().get(),&this->_groupRank);

      assert( this->_groupRank == (globalRank % _physics->comm().size() ) );


      //check against available physics solutions
      std::set<std::string> physicsSolutions = physics->getAvailableSolutions();
      if (this->_solutionNames.empty())
        this->_solutionNames = physicsSolutions ;
      else
      {
        std::set<std::string>::iterator solName = this->_solutionNames.begin() ;
        for( ; solName != this->_solutionNames.end(); solName++)
          if ( !physicsSolutions.count(*solName) ) 
          {
            std::cerr << std::endl;
            std::cerr 
              << "ERROR: requested solution "
              << *solName << " " 
              << " is not available in this PhysicsModel class. \n" ;
            std::cerr << std::endl;
            std::abort();
          }
      }


    }

/********************************************//*
 * \brief Secondary constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
          std::shared_ptr<SurrogateModel<T_S,T_P> > primarySurrogate, 
          std::vector<unsigned int> increaseOrder ,
          unsigned int multiplyOrder ,
          std::set<std::string> evaluateSolutions , 
          std::set<std::string> computeSolutions 
        )
      : 
        _physics( primarySurrogate->getPhysics() ), 
        _evalSurrogate( primarySurrogate ),
        SurrogateModelBase<T_S,T_P>(
            primarySurrogate->getComm(),
            primarySurrogate->getParameters(),
            primarySurrogate->getExpansionOrder(),
            computeSolutions)
    {
      int globalRank,globalSize;
      MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
      MPI_Comm_size(MPI_COMM_WORLD,&globalSize);
      this->_physicsGroup = globalRank / ( _physics->comm().size() ) ;
      this->_nPhysicsGroups = globalSize / ( _physics->comm().size() ) ;
      MPI_Comm_rank(_physics->comm().get(),&this->_groupRank);

      assert( this->_groupRank == (globalRank % _physics->comm().size() ) );


      // augment order appropriately
      this->_order = primarySurrogate->getExpansionOrder() ;
      if (increaseOrder.empty())
        _increaseOrder = std::vector<unsigned int>(this->_dimension,0) ;
      else
        _increaseOrder = increaseOrder;

      agnos_assert(_increaseOrder.size() == this->_dimension);

      // check/reset increase order for CONSTANT parameters
        for( unsigned int i=0; i<this->_parameters.size(); i++ )
          if ( this->_parameters[i]->type() == CONSTANT )
            if (_increaseOrder[i] != 0)
            {
              std::cout 
                << "WARNING: forcing CONSTANT parameter (" 
                << i 
                << ") increaseOrder to 0 \n" ;
              _increaseOrder[i] = 0;
            }

      _multiplyOrder = multiplyOrder ;
      for (unsigned int i=0; i<this->_order.size();i++)
      {
        this->_order[i] += _increaseOrder[i] ;
        this->_order[i] *= _multiplyOrder ;
      }

      // if no output solutions are provided assume we want all 
      if (computeSolutions.size() == 0)
        this->_solutionNames = primarySurrogate->getSolutionNames();
      else
        this->_solutionNames  = computeSolutions ;

      // if no evaluateSolutions are provided assume we use all
      if (evaluateSolutions.size() == 0)
        _evalNames = this->_solutionNames ;
      else
        _evalNames = evaluateSolutions ;
      

    }


/********************************************//**
 * \brief Default destructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::~SurrogateModel( )
    {
    }



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::refine( )
    {
      std::vector<unsigned int> increase(this->_dimension, 1);
      refine( increase );
    }


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::refine( 
        const std::vector<unsigned int>& increase 
        )
    {

      // verify correct dimension of increase is given
      agnos_assert( (increase.size() == this->_parameters.size()) );

      // if this is a primary surrogate incease order
      if ( _evalSurrogate == NULL )
      {
        for( unsigned int i=0; i<this->_parameters.size(); i++ )
        {
          // increase order
          this->_order[i] += increase[i] ;

          // if its a constant parameter and we raised the order reset it back
          // to zero
          if ( this->_parameters[i]->type() == CONSTANT )
            if (this->_order[i] != 0)
            {
              std::cout 
                << "WARNING: forcing CONSTANT parameter (" 
                << i 
                << ") order to 0th order \n" ;
              this->_order[i] = 0;
            } //end if increase != 0
          //end if paramType = CONSTANT
          
        } //end for each parameter
      }
      // otherwise refine based on primary surrogate
      else
      {
        this->_order = _evalSurrogate->getExpansionOrder();
        for (unsigned int i=0; i<this->_order.size();i++)
        {
          this->_order[i] += _increaseOrder[i] ;
          this->_order[i] *= _multiplyOrder ;
        }

      }

      /* this->initialize(); */
      /* this->build(); */
    }

  template class 
    SurrogateModel<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;


}
