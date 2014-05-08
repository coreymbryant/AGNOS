
// local
#include "agnosDefines.h"
#include "EvaluatorPseudoSpectral.h"
#include "PhysicsCatenary.h"
#include "PseudoSpectralTensorProduct.h"

// boost
#define BOOST_TEST_MODULE evaluatorTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace AGNOS;

BOOST_AUTO_TEST_CASE( evaluator_constructor )
{
  
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "dummy");
  
  av[0] = name;
  
  // Initialize libmesh
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  PETSC_COMM_WORLD = MPI_COMM_WORLD ;
  int ierr = PetscInitialize(&ac, const_cast<char***>(&av),NULL,NULL);
  LibMeshInit init (ac, av);

  // construct parameter object
  unsigned int dimension = 3;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int p=0;p<dimension;p++)
    parameters.push_back(
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("UNIFORM",0.0,1.0) )
        );

  // order
  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(0);

  // compute solutions names
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("qoi");

  SurrogateEvaluator<T_S,T_P>* evaluator =
    new EvaluatorPseudoSpectral<T_S,T_P>(
        comm,
        parameters,
        order,
        computeSolutions
        );

  BOOST_REQUIRE( ( evaluator->dimension() == parameters.size() ) ) ;

  std::vector<unsigned int> evalOrder = evaluator->getExpansionOrder();
  BOOST_REQUIRE( ( evalOrder.size() == order.size() ) ) ;
  for(unsigned int i=0; i<order.size();i++)
    BOOST_REQUIRE( ( evalOrder[i] == order[i] ) ) ;

  std::set<std::string> evalSolutions = evaluator->getSolutionNames();
  BOOST_REQUIRE( ( evalSolutions.size() == computeSolutions.size() ) ) ;
  std::set<std::string>::iterator sit = evalSolutions.begin();
  for(;sit!=evalSolutions.end();sit++)
    BOOST_REQUIRE( ( computeSolutions.count(*sit) != 0 ) ) ;

}

BOOST_AUTO_TEST_CASE( evaluator_copy_constructor )
{
  
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "test");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  // order
  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(0);
  order.push_back(0);

  // compute solutions names
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("qoi");

  // construct a physics object
  GetPot inputfile = GetPot() ;
  std::shared_ptr<AGNOS::PhysicsModel<T_S,T_P> > physics(
      new AGNOS::PhysicsCatenary<T_S,T_P>(comm,inputfile ) 
      );
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int i=0; i<order.size();i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter(UNIFORM, 1.0,3.0) )
        );

  // construct a SurrogateModel
  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> > surrogate (
      new AGNOS::PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        physics,
        parameters,
        order,
        computeSolutions
        )
      );
  std::cout << "pre build" << std::endl;
  surrogate->build();
  std::cout << "post build" << std::endl;
  std::map< std::string, LocalMatrix >  coefficients
    = surrogate->getCoefficients();
  std::vector< std::vector<unsigned int> > indexSet
    = surrogate->indexSet();


  // construct an evaluator based on this surrogate
  std::shared_ptr<SurrogateEvaluator<T_S,T_P> > evaluator(
    new EvaluatorPseudoSpectral<T_S,T_P>( *surrogate ) 
    );

  std::map< std::string, LocalMatrix >  evalCoefficients
    = evaluator->getCoefficients();
  std::vector< std::vector<unsigned int> > evalIndexSet
    = evaluator->indexSet();

  BOOST_REQUIRE( ( evalCoefficients.size() == coefficients.size() ) ) ;
  std::map< std::string, LocalMatrix >::iterator  cit 
    = evalCoefficients.begin();
  for(; cit!=evalCoefficients.end();cit++)
  {
    BOOST_REQUIRE( ( coefficients.count(cit->first) != 0 ) ) ;
    BOOST_REQUIRE( 
        ( coefficients[cit->first].m() ==  cit->second.m() ) ) ;
    BOOST_REQUIRE( 
        ( coefficients[cit->first].n() ==  cit->second.n() ) ) ;
    for(unsigned int i=0; i<cit->second.m();i++)
      for(unsigned int j=0; j<cit->second.n();j++)

        BOOST_REQUIRE( ( cit->second(i,j) == coefficients[cit->first](i,j) ) ) ;
  }

  BOOST_REQUIRE( ( evalIndexSet.size() == indexSet.size() ) ) ;
  for(unsigned int i=0; i<indexSet.size();i++)
  {
    BOOST_REQUIRE( ( evalIndexSet[i].size() == indexSet[i].size() ) ) ;
    for(unsigned int j=0; j<indexSet[i].size();j++)
      BOOST_REQUIRE( ( evalIndexSet[i][j] == indexSet[i][j] ) ) ;
  }



}

BOOST_AUTO_TEST_CASE( evaluator_build )
{
  
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "test");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  // order
  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(0);
  order.push_back(0);

  // compute solutions names
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("qoi");

  // construct a physics object
  GetPot inputfile = GetPot() ;
  std::shared_ptr<AGNOS::PhysicsModel<T_S,T_P> > physics(
      new AGNOS::PhysicsCatenary<T_S,T_P>(comm,inputfile ) 
      );
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int i=0; i<order.size();i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter(UNIFORM, 1.0,3.0) )
        );

  // construct a SurrogateModel
  std::shared_ptr<AGNOS::SurrogateModel<T_S,T_P> > surrogate (
      new AGNOS::PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        physics,
        parameters,
        order,
        computeSolutions
        )
      );
  surrogate->build();
  std::map< std::string, LocalMatrix >  coefficients
    = surrogate->getCoefficients();
  std::vector< std::vector<unsigned int> > indexSet
    = surrogate->indexSet();


  // construct an evaluator based on this surrogate
  std::shared_ptr<SurrogateEvaluator<T_S,T_P> > evaluator(
    new EvaluatorPseudoSpectral<T_S,T_P>(
        comm,
        parameters,
        order,
        computeSolutions
        )
    );

  evaluator->build(indexSet,coefficients);
  std::map< std::string, LocalMatrix >  evalCoefficients
    = evaluator->getCoefficients();
  std::vector< std::vector<unsigned int> > evalIndexSet
    = evaluator->indexSet();

  BOOST_REQUIRE( ( evalCoefficients.size() == coefficients.size() ) ) ;
  std::map< std::string, LocalMatrix >::iterator  cit 
    = evalCoefficients.begin();
  for(; cit!=evalCoefficients.end();cit++)
  {
    BOOST_REQUIRE( ( coefficients.count(cit->first) != 0 ) ) ;
    BOOST_REQUIRE( 
        ( coefficients[cit->first].m() ==  cit->second.m() ) ) ;
    BOOST_REQUIRE( 
        ( coefficients[cit->first].n() ==  cit->second.n() ) ) ;
    for(unsigned int i=0; i<cit->second.m();i++)
      for(unsigned int j=0; j<cit->second.n();j++)

        BOOST_REQUIRE( 
            ( cit->second(i,j) == coefficients[cit->first](i,j) ) ) ;
  }

  BOOST_REQUIRE( ( evalIndexSet.size() == indexSet.size() ) ) ;
  for(unsigned int i=0; i<indexSet.size();i++)
  {
    BOOST_REQUIRE( ( evalIndexSet[i].size() == indexSet[i].size() ) ) ;
    for(unsigned int j=0; j<indexSet[i].size();j++)
      BOOST_REQUIRE( ( evalIndexSet[i][j] == indexSet[i][j] ) ) ;
  }



}

