#ifndef PHYSICS_LIBMESH_H
#define PHYSICS_LIBMESH_H

/* #include "agnosDefines.h" */
#include "PhysicsModel.h"

// libmesh includes
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/exodusII_io.h"

#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_UNORDERED_SET

namespace AGNOS
{


  /********************************************//**
   * \brief Base libmesh physics model class. All libmesh physics, and maybe
   * GRINS, PhysicsModel classes should be derived from here.
   *
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsLibmesh : public PhysicsModel<T_S,T_P>
  {

    public:

      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsLibmesh( const Communicator& comm_in, const GetPot& input );

      /** Destructor */
      virtual ~PhysicsLibmesh( );

      /** Function called by SurrogateModel to solve for requested solution
       * vectors at each evaluation point in parameter space  */
      virtual void compute( 
          std::set<std::string>& computeSolutions,
          const T_S& paramVector, 
          std::map<std::string, T_P >& solutionVectors 
          ) ;

      /** Refinement methods for the physics model */
      void refine (  ) ;
      void refine ( T_P& errorIndicators ) ;

      /** Exact qoi */
      virtual T_P exactQoi( ) 
      {
        T_P resultVector ;
        return resultVector ;
      }

      /** Return libMesh system object */
      const libMesh::System& getSystem( ) const { return *_system; }
      /** Return libMesh es object */
      const libMesh::EquationSystems& getEquationSystems( ) const { 
        return *_equationSystems; }
      /** Return libMesh mesh object */
      const libMesh::MeshBase& getMesh( ) const { return *_mesh; }
      /** Return libMesh mesh refinement object */
      const libMesh::MeshRefinement& getMeshRefinement( ) const { 
        return *_meshRefinement; }
      /** Return libMesh estimator object */
      const libMesh::AdjointRefinementEstimator& getEstimator( ) const { 
        return *_errorEstimator; }
      /** Return libMesh QoISet object */
      const libMesh::QoISet& getQois( ) const { return *_qois; }

    protected:
      /** mesh and equation pointers */
      libMesh::System*                      _system;
      libMesh::EquationSystems*             _equationSystems;
      libMesh::MeshBase*                    _mesh; 
      libMesh::MeshRefinement*              _meshRefinement;
      libMesh::AdjointRefinementEstimator*  _errorEstimator;
      libMesh::QoISet*                      _qois;

      /** refinement options */
      bool         _useUniformRefinement;
      bool         _resolveAdjoint;
      unsigned int _numberHRefinements;
      unsigned int _numberPRefinements;
      unsigned int _maxRefineSteps;

      /** control vizualization output */
      bool         _writePrimalViz;
      bool         _writeAdjointViz;
      virtual void _outputPrimalViz( const T_S& paramVector );
      virtual void _outputAdjointViz( const T_S& paramVector );

      /** build mesh refinement object (can be overidden in derived class) */
      virtual void _buildMeshRefinement();

      /** build error estimator object. Defaults to adjoint_residual_estimator
       * (can be overridden in derived class) */
      virtual void _buildErrorEstimator();

      /** Insert solution vector to map. Helper function used in compute(..)*/
      void _insertSolVec( 
          libMesh::NumericVector<Number>& vec,
          std::string vecName,
          std::map<std::string, T_P >& mapOfVecs 
          ) ;
      /** Insert solution vector to map. Helper function used in compute(..)*/
      void _insertSolVec( 
        T_P& vec, std::string vecName, std::map<std::string, T_P >& mapOfVecs ) 
      {
        mapOfVecs.insert( 
            std::pair<std::string,T_P>( vecName, vec )
              );
        return;
      }


      /** Solve primal problem. Called in compute(..) if requested */
      virtual void _solve( ) ;

      /** Solve adjoint problem. Called in compute(..) if requested */
      virtual void _adjointSolve( libMesh::QoISet& qois ) ;

      /** Set primal solution to provided vector */
      virtual void _setPrimalSolution( T_P& solutionVector );

      /** Set adjoint solution 'j' to provided vector */
      virtual void _setAdjointSolution( 
          T_P& solutionVector, unsigned int j=0 );

      /** Reference to qoi value */
      virtual libMesh::Number _getQoIValue( ) { return _system->qoi[0]; }

      /** dummy timestep for printing visualizations based on parameter values */
      unsigned int _timeStep;

      /** derived PhysicsModel classes need to handle settig parameter values
       * themselves */
      virtual void _setParameterValues( const T_S& parameterValues ) 
      {
        std::cout << std::endl ;
        std::cout 
          << "WARNING: _setParameterValues should be implemented by derived\n"
          << "         classes. Without this interface AGNOS can not change\n"
          << "         system parameters and properly construct a surrogate\n"
          << "         model.\n"
          << std::endl ;
        return;
      }

      /** output file objects */
      libMesh::ExodusII_IO* _adjointExio;
      libMesh::ExodusII_IO* _primalExio;

      /** initialize output objects */
      void _initOutput( )
      {
        if ( _mesh != NULL)
        {
          _adjointExio = new libMesh::ExodusII_IO(*_mesh);
          _primalExio = new libMesh::ExodusII_IO(*_mesh);
        }
        return;
      }

      /** output visualization in exodus format */
      void _writeExodus(
        EquationSystems &es,
        const T_S& paramVector,       // The parameter vector
        libMesh::ExodusII_IO* exio,  // primal or adjoint solve
        std::string fileName);

      /** output visualization in gnuplot format */
      void _writeGnuplot(
        EquationSystems &es,
        const T_S& paramVector,       // The parameter vector
        std::string solution_type = "primal"); // primal or adjoint solve


  };




}

#endif // PHYSICS_LIBMESH_H
