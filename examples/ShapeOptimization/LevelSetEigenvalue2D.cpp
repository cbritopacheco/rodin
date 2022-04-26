#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>
#include <vector>
#include <typeinfo>

// Dependencies required for eigenvalue computation
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsShiftSolver.h>

using namespace Rodin;
using namespace Rodin::Variational;
using namespace Rodin::External;

Eigen::SparseMatrix<double> mfemToEigenSparse(mfem::SparseMatrix m)
{
  typedef Eigen::Triplet<double> Triplet;

  int nonZero = m.NumNonZeroElems();
  std::vector<Triplet> tripletList;

  // Get the values (i,j, data)
  int* m_I = m.GetI();
  int* m_J = m.GetJ();
  double* m_data = m.GetData();

  for(int r = 0; r < m.NumRows(); ++r){
    for(int j=m_I[r]; j<m_I[r+1]; ++j){
      tripletList.push_back(Triplet(r, m_J[j], m_data[j])); // ???? Inverse row and column since Eigen is CSC
    }
  }
  Eigen::SparseMatrix<double> m_sparse(m.NumRows(), m.NumCols());
  m_sparse.setFromTriplets(tripletList.begin(), tripletList.end()); //SEGFAULT

  return m_sparse;
}


class EigenSolver{

  public:
    // Constructors
    EigenSolver();

    // Destructor
    ~EigenSolver() = default;

    // Accessors
    GridFunction<H1>& getEigenFunction(int index);

    double getEigenValue(int index);
    std::vector<double> getEigenValues();

    // Options
    EigenSolver& setNumEV(unsigned int nev);
    EigenSolver& setShift(double shift);

    // Solves the problem
    void solve(mfem::SparseMatrix A, mfem::SparseMatrix B, FiniteElementSpace<H1>& fes);

  private:
    std::vector<double> m_eigenvalues;
    std::vector<GridFunction<H1>*> m_eigenfunctions;
    unsigned int m_nev;
    double m_shift;
};



int main(int argc, char** argv)
{
  const char* meshFile = "../resources/mfem/levelset-cantilever2d-example.mesh";

  // Define interior and exterior for level set discretization
  int Interior = 1, Exterior = 2;

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  // Solver for hilbertian regularization
  auto solver = Solver::UMFPack();

  // Optimization parameters
  size_t maxIt = 10;
  double eps = 1e-6;
  double hmax = 0.05;
  int k = 1; // The eigenvalue to optimize
  auto alpha = ScalarFunction(4 * hmax * hmax); // Parameter for hilbertian regularization

  std::vector<double> obj;

  /*
  // Scalar field finite element space over the whole domain
  FiniteElementSpace<H1> Vh(Omega);

  // Trim the exterior part of the mesh to solve the elasticity system
  SubMesh trimmed = Omega.trim(Exterior, Gamma);

  // Build a finite element space over the trimmed mesh
  FiniteElementSpace<H1> VhInt(trimmed);

  // Elasticity equation
  TrialFunction uInt(VhInt);
  TestFunction  vInt(VhInt);

  // A
  Problem stiffness(uInt, vInt);
  stiffness = Integral(Grad(uInt), Grad(vInt));
  stiffness.update().assemble();

  // B
  Problem mass(uInt, vInt);
  mass = Integral(uInt, vInt);
  mass.update().assemble();

  auto& m1 = stiffness.getStiffnessMatrix();
  auto& m2 = mass.getStiffnessMatrix();


  // Solve eigenvalue problem
  EigenSolver ES;
  ES.setShift(1.0).setNumEV(k+4).solve(m1, m2, VhInt);

  // Get solution and transfer to original domain
  GridFunction u(Vh);
  ES.getEigenFunction(k).transfer(u);
  double mu = ES.getEigenValue(k);

  // Save solution
  u.save("u.gf");
  Omega.save("Omega.mesh");
  */

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    // Scalar field finite element space over the whole domain
    FiniteElementSpace<H1> Vh(Omega);

    // Trim the exterior part of the mesh to solve the elasticity system
    SubMesh trimmed = Omega.trim(Exterior, Gamma);

    // Build a finite element space over the trimmed mesh
    FiniteElementSpace<H1> VhInt(trimmed);

    // Elasticity equation
    TrialFunction uInt(VhInt);
    TestFunction  vInt(VhInt);

    // A
    Problem stiffness(uInt, vInt);
    stiffness = Integral(Grad(uInt), Grad(vInt));
    stiffness.update().assemble();

    // B
    Problem mass(uInt, vInt);
    mass = Integral(uInt, vInt);
    mass.update().assemble();

    auto& m1 = stiffness.getStiffnessMatrix();
    auto& m2 = mass.getStiffnessMatrix();


    // Solve eigenvalue problem
    EigenSolver ES;
    ES.setShift(1.0).setNumEV(k+4).solve(m1, m2, VhInt);

    // Get solution and transfer to original domain
    GridFunction u(Vh);
    ES.getEigenFunction(k).transfer(u);
    double mu = ES.getEigenValue(k);

    // Compute the boundary shape gradientv
    // Vector field finite element space over the whole domain
    FiniteElementSpace<H1> Uh(Omega, 2);
    auto n = Normal(2);
    GridFunction dJ(Uh);
    //dJ = Omega.getVolume(Interior)*(Dot(Grad(u), Grad(u)) + mu*Dot(u,u))*n + mu*n; //  /Integral(Dot(u,u)) sur le sous-mesh
    dJ = Dot(Grad(u), Grad(u))*n;
    //std::cout << dJ.getFiniteElementSpace().getVectorDimension() << std::endl;


    // Hilbert extension-regularization procedure
    // This is bullshit
    TrialFunction g(Uh);
    TestFunction  v(Uh);
    Problem hilbert(g, v);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
            + Integral(g, v)
            + DirichletBC(g, VectorFunction(dJ)).on(GammaN);
    solver.solve(hilbert);

    // Update objective
    obj.push_back(
        mu * Omega.getVolume(Interior));
    std::cout << "[" << i << "] Objective: " << obj.back() << std::endl;

    // Convert data types to mmg types
    auto mmgMesh = Cast(Omega).to<MMG::Mesh2D>();
    auto mmgVel = Cast(g.getGridFunction()).to<MMG::VectorSolution2D>(mmgMesh);

    // Generate signed distance function
    auto mmgLs = MMG::Distancer2D().setInteriorDomain(Interior).distance(mmgMesh);

    // Advect the level set function
    double gInf = std::max(g.getGridFunction().max(), -g.getGridFunction().min());
    double dt = 4 * hmax / gInf;
    MMG::Advect2D(mmgLs, mmgVel).step(dt);

    // Recover the implicit domain
    auto mmgImplicit =
      MMG::ImplicitDomainMesher2D().split(Interior, {Interior, Exterior})
                                   .split(Exterior, {Interior, Exterior})
                                   .setRMC(1e-3)
                                   .setHMax(hmax)
                                   .setBoundaryReference(Gamma)
                                   .discretize(mmgLs);

    // Convert back to Rodin data type
    Omega = Cast(mmgImplicit).to<Rodin::Mesh<>>();

    // Save mesh
    Omega.save("Omega.mesh");


    // Test for convergence
    if (obj.size() >= 2 && abs(obj[i] - obj[i - 1]) < eps)
    {
      std::cout << "Convergence!" << std::endl;
      break;
    }

    std::ofstream plt("obj.txt", std::ios::trunc);
    for (size_t i = 0; i < obj.size(); i++)
      plt << i << "," << obj[i] << "\n";
  }



  return 0;
}


// Constructor
EigenSolver::EigenSolver(){
  m_nev = 0;
  m_shift = 0.0;
}

// Accessors
GridFunction<H1>& EigenSolver::getEigenFunction(int index){
  return *m_eigenfunctions[index];
}

double EigenSolver::getEigenValue(int index){
    return m_eigenvalues[index];
}

std::vector<double> EigenSolver::getEigenValues(){
  return m_eigenvalues;
}

// Options
EigenSolver& EigenSolver::setNumEV(unsigned int nev){
  m_nev = nev;
  return *this;
}

EigenSolver& EigenSolver::setShift(double shift){
  m_shift = shift;

  return *this;
}

void EigenSolver::solve(mfem::SparseMatrix A,mfem::SparseMatrix B, FiniteElementSpace<H1>& fes){

  // 1. CA and B to Eigen::SparseMatrix
  //WARNING: MAY BE USELESS SINCE SPECTRA HANDLE OTHER MATRICES TYPES
  //BUT MAY BE FASTER. IN FACT I DON'T HAVE A CLUE.
  Eigen::SparseMatrix<double> A_sparse = mfemToEigenSparse(A);
  Eigen::SparseMatrix<double> B_sparse = mfemToEigenSparse(B);

  // 2. Use Spectra with the shift-inverse method to compute eigenvalue
  // Define the shift-inverse and vector-multiplication operations
  using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  using BOpType = Spectra::SparseSymMatProd<double>;
  OpType op(A_sparse, B_sparse);
  BOpType Bop(B_sparse);

  // Construct the generalized eigensolver object
  int ncv = 2*m_nev+1;
  Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
    geigs(op, Bop, m_nev, ncv, m_shift);


  // Initialize and compute
  geigs.init();
  int nconv = geigs.compute(Spectra::SortRule::LargestMagn); // UTILISER LARGESTMAGN ET PAS SMALLESTMAGN CAR ON EST EN SHIFT-INVERSE MODE

  // Retrieve results
  Eigen::VectorXd evalues;
  Eigen::MatrixXd evecs;
  if (geigs.info() == Spectra::CompInfo::Successful)
  {
      evalues = geigs.eigenvalues();
      evecs = geigs.eigenvectors();
  }


  // Store in m_eigenvalues and m_eigenfunctions
  // (and reverse the order to order from lowest to higest)
  int dim = evecs.rows();

  for(int i = 0; i < m_nev; i++){
    m_eigenvalues.push_back(evalues[m_nev-i-1]);

    std::unique_ptr<double[]> data(new double[dim]);
    for(int j=0; j<dim; j++){
      data.get()[j] = evecs.col(m_nev-i-1).data()[j];
    }
    //std::copy(evecs.col(m_nev-i-1).data(), evecs.col(m_nev-i-1).data()+dim, data.get());

    m_eigenfunctions.push_back(new GridFunction<H1>(fes));
    m_eigenfunctions[i]->setData(move(data), dim);
  }
}
