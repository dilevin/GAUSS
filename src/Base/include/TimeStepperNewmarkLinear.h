#ifndef TimeStepperNewmarkLinear_h
#define TimeStepperNewmarkLinear_h

#include <World.h>
#include <Assembler.h>
#include <TimeStepper.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <UtilitiesEigen.h>
#include <UtilitiesMATLAB.h>
#include <Eigen/SparseCholesky>
#include <SolverPardiso.h>

namespace Gauss {

	template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
	class TimeStepperImpNewmarkLinear {
	public:

		TimeStepperImpNewmarkLinear(Eigen::SparseMatrix<double> &P) {
			m_P = P;
			m_factored = false;
		}

		TimeStepperImpNewmarkLinear(const TimeStepperImpNewmarkLinear &toCopy) {
			m_factored = false;
		}

		~TimeStepperImpNewmarkLinear() {
		}

		//Methods
		template<typename World>
		void step(World &world, double dt, double t);

		inline typename VectorAssembler::MatrixType & getLagrangeMultipliers() { return m_lagrangeMultipliers; }

	protected:

		bool initialized = false;

		MatrixAssembler m_massMatrix;
		MatrixAssembler m_stiffnessMatrix;

		Eigen::SparseMatrix<double> m_P;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

		//storage for lagrange multipliers
		typename VectorAssembler::MatrixType m_lagrangeMultipliers;

		bool m_factored, m_refactor;
	};
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
template<typename World>
void TimeStepperImpNewmarkLinear<DataType, MatrixAssembler, VectorAssembler>::step(World &world, double dt, double t) {
	// precompute and prefactor Hessian
	if (!initialized) {
		getMassMatrix(m_massMatrix, world);
		getStiffnessMatrix(m_stiffnessMatrix, world);

		Eigen::SparseMatrix<DataType, Eigen::RowMajor> H;
		H = m_P*(*m_massMatrix)*m_P.transpose() - (dt*dt/4.0)*m_P*(*m_stiffnessMatrix)*m_P.transpose();
		solver.compute(H);
		initialized = true;
	}
	// assemble b vector
	Eigen::Map<Eigen::VectorXd> q = mapStateEigen<0>(world);
	Eigen::Map<Eigen::VectorXd> qDot = mapStateEigen<1>(world);
	AssemblerEigenVector<double> force;
	getForceVector(force, world);
	Eigen::VectorXd  b = dt * m_P*(*m_massMatrix)*m_P.transpose()*m_P*qDot + (dt*dt / 2.0)*m_P*(*force);
	// solve for delta
	Eigen::VectorXd delta = m_P.transpose() * solver.solve(b);
	// update state
	q = q + delta;
	qDot = (2.0 / dt)*delta - qDot;
}

template<typename DataType, typename MatrixAssembler, typename VectorAssembler>
using TimeStepperNewmarkLinear = TimeStepper<DataType, TimeStepperImpNewmarkLinear<DataType, MatrixAssembler, VectorAssembler> >;

#endif //TimeStepperNewmarkLinear_h