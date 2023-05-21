#ifndef CORENC_SOLVERS_DG_SOLVER_H_
#define CORENC_SOLVERS_DG_SOLVER_H_

#include "../CoreNCFEM/Grids/TriangularMesh.h"
#include "../CoreNCFEM/Methods/DGMethod.h"
#include "../Problems/DiffusionScalar.h"
#include "../CoreNCA/MatrixSkyline.h"
#include "../CoreNCFEM/Methods/FEAnalysis.h"

namespace corenc
{
	namespace solvers
	{
		template<class _Problem, class _Mesh, class _Result>
        class dg_solver
		{
            using _Method = method::DGMethod<_Problem, _Mesh, Algebra::MatrixSkyline>;
		public:
            dg_solver() :m_method{ nullptr } {}
            ~dg_solver()
			{
				if (m_method != nullptr)
					delete m_method;
			}
			// terms, method, mesh, solver, result
			const int				elliptic_solver(_Problem*, _Mesh*, _Result*);
            const double			get_value(const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const double			get_value(const _Method*, const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const double			get_value(const _Mesh&, const _Result&, const Mesh::Point& p, const int i) const;
            const Mesh::Point		get_gradvalue(const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const Mesh::Point		get_gradvalue(const _Mesh&, const _Result&, const Mesh::Point& p, const int i) const;
		private:
			_Method * m_method;
		};

		template<class _Problem, class _Mesh, class _Result>
        const int dg_solver<_Problem, _Mesh, _Result>::elliptic_solver(_Problem* problem, _Mesh* mesh, _Result* result)
		{
			std::vector<double> res;
			//std::shared_ptr<Algebra::MatrixSkyline> matrix{ new Algebra::MatrixSkyline() };
			Algebra::MatrixSkyline* matrix{ new Algebra::MatrixSkyline() };
			std::vector<double> rhs;
			if (m_method != nullptr)
				delete m_method;
			m_method = new _Method{ problem, mesh, matrix, &rhs };
			m_method->Discretization();
			Algebra::ESolver esl{ Algebra::Solvers::GMRES };
			*result = esl.Solve(*matrix, rhs, *result, res, 100000, 1e-13);
			delete matrix;
			return 0;
		}

		template<class _Problem, class _Mesh, class _Result>
        const double dg_solver<_Problem, _Mesh, _Result>::get_value(const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
		{
			if (m_method != nullptr)
				return m_method->GetSolution(mesh, res, p);
			return 0.;
		}

		template<class _Problem, class _Mesh, class _Result>
        const Mesh::Point dg_solver<_Problem, _Mesh, _Result>::get_gradvalue(const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
		{
			if (m_method != nullptr)
				return m_method->GetGradSolution(mesh, res, p);
			return Mesh::Point(0, 0, 0);
		}

        template<class _Problem, class _Mesh, class _Result>
        const double dg_solver<_Problem, _Mesh, _Result>::get_value(const _Method* method2, const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
        {
            if (method2 != nullptr)
                return method2->GetSolution(mesh, res, p);
            return 0.;
        }

        template<class _Problem, class _Mesh, class _Result>
        const double dg_solver<_Problem, _Mesh, _Result>::get_value(const _Mesh& mesh, const _Result& res, const Mesh::Point& p, const int i) const
        {
            if (m_method != nullptr)
                return m_method->GetSolution(mesh, res, p, i);
            return 0.;
        }
        template<class _Problem, class _Mesh, class _Result>
        const Mesh::Point dg_solver<_Problem, _Mesh, _Result>::get_gradvalue(const _Mesh& mesh, const _Result& res, const Mesh::Point& p, const int i) const
        {
            if (m_method != nullptr)
                return m_method->GetGradSolution(mesh, res, p, i);
            return Mesh::Point(0, 0, 0);
        }

	}
}
#endif // !CORENC_SOLVERS_dg_solver_H_
