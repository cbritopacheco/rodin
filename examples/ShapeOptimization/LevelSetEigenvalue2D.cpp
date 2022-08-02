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

static const char* meshFile = "../resources/mfem/levelset-cantilever2d-example.mesh";

// Define interior and exterior for level set discretization
static constexpr int Interior = 1, Exterior = 2;

// Define boundary attributes
static constexpr int Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Optimization parameters
static constexpr size_t maxIt = 36;
static constexpr double eps = 1e-6;
static constexpr double hmax = 0.05;
static constexpr double target_volume = 0.8;
static constexpr int k = 2; // The eigenvalue to optimize
static constexpr double alpha = 4 * hmax * hmax; // Parameter for hilbertian regularization
static constexpr double ell = 1.0;

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


class EigenSolver
{

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
  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  // Solver for hilbertian regularization
  auto solver = Solver::UMFPack();

  std::ofstream plt("obj.txt", std::ios::trunc);

  // Optimization loop
  std::vector<double> obj;
  for (size_t i = 0; i < maxIt; i++)
  {
    // Scalar field finite element space over the whole domain
    FiniteElementSpace<H1> Vh(Omega);

    // Trim the exterior part of the mesh to solve the elasticity system
    SubMesh trimmed = Omega.trim(Exterior);

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

    // Hilbert extension-regularization procedure
    auto n = Normal(2);

    auto gu = Grad(u);
    gu.traceOf(Interior);

    auto dJ = Dot(gu, gu) * n
            - mu * u * u * n
            - 2 * ell * (Omega.getVolume(Interior) - target_volume) * n;

    FiniteElementSpace<H1> Uh(Omega, 2);
    TrialFunction g(Uh);
    TestFunction  v(Uh);
    Problem hilbert(g, v);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
            + Integral(g, v)
            - BoundaryIntegral(dJ, v).over(Gamma);
    solver.solve(hilbert);

    // Update objective
    obj.push_back(
        mu * Omega.getVolume(Interior));
    std::cout << "[" << i << "] Objective: " << obj.back() << std::endl;

    // Generate signed distance function
    FiniteElementSpace<H1> Dh(Omega);
    auto dist = MMG::Distancer(Dh).setInteriorDomain(Interior)
                                  .distance(Omega);

    // Advect the level set function
    double gInf = std::max(g.getGridFunction().max(), -g.getGridFunction().min());
    double dt = 4 * hmax / gInf;
    MMG::Advect(dist, g.getGridFunction()).step(dt);

    // Recover the implicit domain
    Omega = MMG::ImplicitDomainMesher().split(Interior, {Interior, Exterior})
                                       .split(Exterior, {Interior, Exterior})
                                       .setRMC(1e-3)
                                       .setHMax(hmax)
                                       .setBoundaryReference(Gamma)
                                       .discretize(dist);

    MMG::MeshOptimizer().setHMax(hmax).optimize(Omega);

    Omega.save("Omega.mesh");

    // Test for convergence
    if (obj.size() >= 2 && abs(obj[i] - obj[i - 1]) < eps)
    {
      std::cout << "Convergence!" << std::endl;
      break;
    }
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
      data[j] = evecs.col(m_nev-i-1)[j];
    }
    //std::copy(evecs.col(m_nev-i-1).data(), evecs.col(m_nev-i-1).data()+dim, data.get());

    m_eigenfunctions.push_back(new GridFunction<H1>(fes));
    m_eigenfunctions[i]->setData(std::move(data), dim);
  }
}
