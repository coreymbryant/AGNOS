
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

#include "agnosDefines.h"
#include "Parameter.h"
#include "PhysicsFunction.h"

namespace AGNOS
{

  enum SurrogateModelType{
    PSEUDO_SPECTRAL_TENSOR_PRODUCT=0,
    PSEUDO_SPECTRAL_SPARSE_GRID,
    PSEUDO_SPECTRAL_MONTE_CARLO,
    COLLOCATION };

  /********************************************//**
   * \brief Base surrogate model class
   *
   * Allows for derivation of Collocation and Pseudospectral surrogate models.
   * Function to construct SurrogateModel for must me defined by providing a
   * PhysicsFunction object. 
   ***********************************************/
  template<class T_S, class T_P>
  class SurrogateModel
  {


    public: 

      // single physics function constructors
      SurrogateModel(
          const Communicator*               comm,
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const unsigned int                order 
          );
      SurrogateModel(
          const Communicator*               comm,
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      // multiple physics function constructors
      SurrogateModel(
          const Communicator*               comm,
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const unsigned int                                  order 
          );
      SurrogateModel(
          const Communicator*               comm,
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const std::vector<unsigned int>&                    order
          );

      SurrogateModel( );           /**< Default constructor */
      virtual ~SurrogateModel( );  /**< Default destructor */

      // surrogate construction and evaluation
      virtual void build( ) = 0; 

      virtual std::map<std::string, T_P> evaluate( 
          std::vector<std::string> solutionNames,  ///< solution to return
          T_S& parameterValues /**< parameter values to evaluate*/
          ) = 0;
      T_P evaluate( 
          std::string solutionName,  ///< solution to return
          T_S& parameterValues     ///< parameter values to evaluate*/
          ) ;
      std::map<std::string, T_P> evaluate( 
          T_S& parameterValues     ///< parameter values to evaluate*/
          ) ;

      virtual void refine( ) = 0;

      std::map< std::string, T_P > mean( ) ;


      // Manipulators
      void setParameters( std::vector<Parameter*> parameters );
      std::vector<Parameter*>   getParameters( ) const;

      void printCoefficients( 
        std::vector<std::string> solutionNames,
        std::ostream& out ) ;
      void printCoefficients( std::string solutionName, std::ostream& out ) ;
      void printCoefficients( std::ostream& out ) ;
      const std::map< std::string, std::vector<T_P> >   getCoefficients() const;

      virtual void printIntegrationWeights( std::ostream& out ) const = 0;
      virtual void printIntegrationPoints( std::ostream& out ) const = 0;
      virtual void printIndexSet( std::ostream& out ) const = 0;
      virtual void prettyPrintIntegrationWeights( ) const = 0;
      virtual void prettyPrintIntegrationPoints( ) const = 0;
      virtual void prettyPrintIndexSet( ) const = 0;

      std::vector<unsigned int> getExpansionOrder( ) const;

    protected: 
      const Communicator* m_comm;
      
      std::vector<unsigned int>                           m_order;  
      std::map< std::string, std::vector<T_P> >           m_coefficients;
      unsigned int                                        m_totalNCoeff;

      std::map< std::string, unsigned int>                m_solSize;

      std::map< std::string, PhysicsFunction<T_S,T_P>* >  m_solutionFunction;
      std::vector<Parameter*>                             m_parameters;
      unsigned int                                        m_dimension;

      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator*               comm,
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        std::vector<Parameter*>                             parameters,
        unsigned int                                        order
        )
      : 
        m_comm(comm),
        m_solutionFunction(solutionFunction), 
        m_parameters(parameters),
        m_dimension( parameters.size() )
    {
      m_order = std::vector<unsigned int>(m_dimension,order);
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator*               comm,
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        std::vector<Parameter*>                   parameters,
        const std::vector<unsigned int>&          order
        )
      : 
        m_comm(comm),
        m_solutionFunction(solutionFunction), 
        m_parameters(parameters),
        m_dimension( parameters.size() ), 
        m_order(order)
    {

      if (m_order.size() != parameters.size() )
      {
        std::cout 
          << std::endl
          << "\tERROR:"
          << " order vector dimension does not match number of parameters"
          << std::endl
          << std::endl;
        assert(0);
      }
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator*               comm,
        PhysicsFunction<T_S,T_P>*                 solutionFunction,
        std::vector<Parameter*>                   parameters,
        const std::vector<unsigned int>&          order
        ) 
      : 
        m_comm(comm),
        m_parameters(parameters), 
        m_dimension( parameters.size() ),
        m_order(order)
    {
      m_solutionFunction.insert( 
          std::pair<std::string, PhysicsFunction<T_S,T_P>* >(
            solutionFunction->name(),solutionFunction) );


      if (m_order.size() != parameters.size() )
      {
        std::cout 
          << std::endl
          << "\tERROR:"
          << " order vector dimension does not match number of parameters"
          << std::endl
          << std::endl;
        assert(0);
      }
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator*               comm,
        PhysicsFunction<T_S,T_P>*                 solutionFunction,
        std::vector<Parameter*>                   parameters,
        unsigned int                              order
        ) 
      : 
        m_comm(comm),
        m_parameters(parameters), 
        m_dimension( parameters.size() )
    {
      m_solutionFunction = std::map< std::string, PhysicsFunction<T_S,T_P>* > ( 
          std::map<std::string, PhysicsFunction<T_S,T_P>* >(
            solutionFunction->name(),solutionFunction) );

      m_order = std::vector<unsigned int>(m_dimension,order);
    }

/********************************************//*
 * \brief Default constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( )
    {
    }

/********************************************//**
 * \brief Default destructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::~SurrogateModel( )
    {
      /* typename std::map< std::string, PhysicsFunction<T_S,T_P>* >::iterator id; */
      /* for (id=m_solutionFunction.begin(); id !=m_solutionFunction.end(); id++) */
      /*   delete id->second; */

      /* for (unsigned int i=0; i<m_parameters.size(); i++) */
      /*   delete m_parameters[i]; */

    }

/********************************************//**
 * \brief set parameters
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::setParameters( 
        std::vector<Parameter*> parameters)
    {
      m_parameters = parameters ;
      m_dimension = parameters.size();
      return;
    }

/********************************************//**
 * \brief get number of parameters
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<Parameter*> SurrogateModel<T_S,T_P>::getParameters( ) const
    {
      return m_parameters ;
    }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<unsigned int> SurrogateModel<T_S,T_P>::getExpansionOrder( )
    const
  {
    return m_order;
  }

/********************************************//**
 * \brief return coefficient values
 ***********************************************/
  template<class T_S, class T_P> 
    const std::map< std::string, std::vector<T_P> >
    SurrogateModel<T_S,T_P>::getCoefficients( ) const
    {
      return m_coefficients;
    }

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::printCoefficients( 
        std::vector<std::string> solutionNames,
        std::ostream& out ) 
    {
      unsigned int myRank = m_comm->rank();
      unsigned int commSize = m_comm->size();
      unsigned int coeffStart = std::min(m_comm->rank(), m_totalNCoeff-1);

      if ( myRank == 0)
        out << "#" << std::string(75,'=') << std::endl;

      for (unsigned int i=0; i < solutionNames.size(); i++)
      {
        std::string id = solutionNames[i];

        if (myRank == 0)
        {
          out << "#" << std::string(75,'-') << std::endl;
          out << "#" << "\t Solution: " << id << std::endl;
          out << "#" << std::string(75,'-') << std::endl;
        }

        for(unsigned int c=0; c<this->m_totalNCoeff; c++)
        {
          unsigned int myC = (c-coeffStart)/commSize;
          std::vector<double> myCoeff(m_solSize[id],0.);
          libMesh::Parallel::MessageTag tag(c);
          libMesh::Parallel::Status stat;

          if (c%commSize == myRank)
          {
            if (myRank == 0)
            {
              for(unsigned int comp=0; comp < (m_coefficients[id])[myC].size(); comp++)
                out << std::setprecision(5) << std::scientific 
                  << (m_coefficients[id])[myC](comp) << " " ;
              out << std::endl;
            }
            else
            {
              for(unsigned int comp=0; comp<myCoeff.size(); comp++)
                myCoeff[comp]= m_coefficients[id][myC](comp) ; 
              m_comm->send(0,myCoeff,tag);

            }
          }
          else
          {
            if (myRank == 0)
            {
              stat=m_comm->receive( c%commSize, myCoeff, tag);
              for(unsigned int comp=0; comp < myCoeff.size(); comp++)
                out << std::setprecision(5) << std::scientific 
                  << myCoeff[comp] << " " ;
              out << std::endl;
            }
          }


        }
      }


    }

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::printCoefficients( 
        std::string solutionName,
        std::ostream& out ) 
    {
      std::vector<std::string> solutionsToPrint(1,solutionName);
      printCoefficients(solutionsToPrint, out);
    }

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::printCoefficients( 
        std::ostream& out ) 
    {
      std::vector< std::string > solutionsToPrint ;
      
      typename std::map< std::string, std::vector<T_P> >::iterator id;
      for (id=m_coefficients.begin(); id!=m_coefficients.end(); id++)
        solutionsToPrint.push_back( id->first ) ;
      printCoefficients(solutionsToPrint, out);
    }


/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S, class T_P> 
    T_P SurrogateModel<T_S,T_P>::evaluate( 
        std::string solutionName,  ///< solution to return
        T_S& parameterValues     ///< parameter values to evaluate*/
        ) 
    {
      std::vector< std::string > solutionsToGet(1,solutionName);
      std::map< std::string, T_P > solutionVectors
        = this->evaluate( solutionsToGet, parameterValues ) ;

      return solutionVectors[solutionName];
    }

/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S, class T_P> 
    std::map<std::string, T_P> SurrogateModel<T_S,T_P>::evaluate( 
        T_S& parameterValues     ///< parameter values to evaluate*/
        ) 
    {
      std::vector< std::string > solutionsToGet ;
      
      typename std::map< std::string, PhysicsFunction<T_S,T_P>* >::iterator id;
      for (id=m_solutionFunction.begin(); id!=m_solutionFunction.end(); id++)
        solutionsToGet.push_back( id->first ) ;

      return evaluate( solutionsToGet, parameterValues ) ;    
    }
  
  /********************************************//**
   * \brief get mean of surrogate model
   *      returns first coeff vector
   * 
   ***********************************************/
  template<class T_S,class T_P>
     std::map< std::string, T_P > SurrogateModel<T_S,T_P>::mean( )
    {
      std::map< std::string, T_P > meanCoefficients;

      typename std::map< std::string, std::vector<T_P> >::iterator id;
      for (id=m_coefficients.begin(); id!=m_coefficients.end(); id++)
        meanCoefficients.insert(
            std::pair< std::string, T_P>( id->first, (id->second)[0])
            );

      return meanCoefficients;
    }



}
#endif //SURROGATE_MODEL_H
