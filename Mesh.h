#ifndef MESH_H
#define MESH_H

#include "Types.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

// Eigen Iterative
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// AMGCL Headers
#include <iostream>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <unsupported/Eigen/SparseExtra> // For reading MatrixMarket files

#include <amgcl/backend/eigen.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>

class Mesh {
public:
    // default constructor
    Mesh();

    // read mesh from file
    bool read(const std::string& fileName);

    // write mesh to file
    bool write(const std::string& fileName) const;

    // computes geodesic distances from source point
    std::vector<double> computeGeodesics(const int vIdx);

    // AMGCL Solver for Geodesics 
    std::vector<double> AMGCLSolver(const int vIdx);

    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Eigen::Vector3d> uvs;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;

private:
    // builds Laplacian operator
    void buildLaplacian(Eigen::SparseMatrix<double , Eigen::RowMajor>& L) const;

    // builds area matrix
    void buildAreaMatrix(Eigen::SparseMatrix<double ,Eigen::RowMajor>& A) const;

    // computes time step
    double computeTimeStep() const;

    // computes gradient for faces
    void computeFaceGradients(Eigen::MatrixXd& gradients, const Eigen::VectorXd& u) const;

    // computes integrated gradient for vertices
    void computeIntegratedDivergence(Eigen::VectorXd& integratedDivs,const Eigen::MatrixXd& gradients) const;

    // center mesh about origin and rescale to unit radius
    void normalize();

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> heatSolver;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> poissonSolver;


};

#endif