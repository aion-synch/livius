#ifndef CORENC_SOLVERS_FEM_SOLVER_H_
#define CORENC_SOLVERS_FEM_SOLVER_H_

#include "../CoreNCFEM/Grids/TriangularMesh.h"
#include "../CoreNCFEM/Methods/FEMethod.h"
#include "../Problems/DiffusionScalar.h"
#include "../CoreNCA/MatrixSkyline.h"
#include "../CoreNCFEM/Methods/FEAnalysis.h"

// FINITE ELEMENT METHOD SOLVER ONLY IN SPACE

namespace corenc
{
	namespace solvers
	{
		template<class _Problem, class _Mesh, class _Result>
		class fem_solver
		{
			using _Method = method::FEMethod<_Problem, _Mesh, Algebra::MatrixSkyline>;
			using _Method2 = method::FEMethod<_Problem, _Mesh, Algebra::Matrix>;
		public:
            fem_solver() :m_method2{ nullptr }, m_method{nullptr}{}
			~fem_solver()
			{
				if(m_method2 != nullptr)
					delete m_method2;
                if(m_method != nullptr)
                    delete m_method;
			}
			// terms, method, mesh, solver, result
            const int				elliptic_solver(_Problem*, _Mesh*, _Result*);
            const int				elliptic_solver_gauss(_Problem*, _Mesh*, _Result*);
			const double			get_value(const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const double			get_value(const _Method2*, const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const double			get_value(const _Method*, const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const double			get_value(const _Mesh&, const _Result&, const Mesh::Point& p, const int i) const;
			const Mesh::Point		get_gradvalue(const _Mesh&, const _Result&, const Mesh::Point& p) const;
            const Mesh::Point		get_gradvalue(const _Mesh&, const _Result&, const Mesh::Point& p, const int i) const;
		private:
            _Method*				m_method;
            _Method*				m_method2;
            //_Method2*				m_method2;
		};

        template<class _Problem, class _Mesh, class _Result>
        const int fem_solver<_Problem, _Mesh, _Result>::elliptic_solver(_Problem* problem, _Mesh* mesh, _Result* result)
        {
            std::vector<double> res;
            std::vector<double> res2;
            //std::shared_ptr<Algebra::MatrixSkyline> matrix{ new Algebra::MatrixSkyline() };
            Algebra::MatrixSkyline* matrix{ new Algebra::MatrixSkyline() };
            std::vector<double> rhs;
            if (m_method != nullptr)
                delete m_method;

            m_method = new _Method{ problem, mesh, matrix, &rhs };

            m_method->Discretization();

            Algebra::ESolver esl{ Algebra::Solvers::BiCGStab };
            //std::cout << "Size:\t" << matrix->GetSize() << std::endl;

            //std::cout << matrix->GetSize() << std::endl;
            *result = esl.Solve(*matrix, rhs, *result, res, 100000, 1e-13);
            //std::cout << matrix->GetSize() << std::endl;

            //esl.Pardiso(*matrix, rhs, *result);
            //res.resize(matrix2->GetSize());
            //for (int i = 0; i < matrix2->GetSize(); ++i)
            //{
                //for (int j = 0; j < matrix2->GetSize(); ++j)
                //{
                    //res[i] += result->operator[](j) * (matrix2->GetElement(i, j));
                //}
            //}
            delete matrix;

            return 0;
        }

		template<class _Problem, class _Mesh, class _Result>
        const int fem_solver<_Problem, _Mesh, _Result>::elliptic_solver_gauss(_Problem* problem, _Mesh* mesh, _Result* result)
		{
			std::vector<double> res;
			std::vector<double> res2;
			//std::shared_ptr<Algebra::MatrixSkyline> matrix{ new Algebra::MatrixSkyline() };
            //Algebra::MatrixSkyline* matrix{ new Algebra::MatrixSkyline() };
            Algebra::Matrix* matrix2{ new Algebra::Matrix() };
			std::vector<double> rhs;
            //if (m_method != nullptr)
            //    delete m_method;
            if (m_method2 != nullptr)
                delete m_method2;
            
            //m_method = new _Method{ problem, mesh, matrix, &rhs };
            m_method2 = new _Method2{ problem, mesh, matrix2, &rhs };
            //m_method->Discretization();
            m_method2->Discretization();
            //Algebra::ESolver esl{ Algebra::Solvers::BiCGStab };
            //std::cout << "Size:\t" << matrix->GetSize() << std::endl;
            Algebra::ESolver esl{ Algebra::Solvers::Gauss };
            //std::cout << matrix->GetSize() << std::endl;
            //*result = esl.Solve(*matrix, rhs, *result, res, 100000, 1e-13);
            //std::cout << matrix->GetSize() << std::endl;
            esl.Gauss(*matrix2, rhs, *result);
            //esl.Pardiso(*matrix, rhs, *result);
			//res.resize(matrix2->GetSize());
			//for (int i = 0; i < matrix2->GetSize(); ++i)
			//{
				//for (int j = 0; j < matrix2->GetSize(); ++j)
				//{
					//res[i] += result->operator[](j) * (matrix2->GetElement(i, j));
				//}
			//}
            //delete matrix;
            delete matrix2;
			return 0;
        }
		template<class _Problem, class _Mesh, class _Result>
		const double fem_solver<_Problem, _Mesh, _Result>::get_value(const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
		{
            if (m_method2 != nullptr)
                return m_method2->GetSolution(mesh, res, p);
			return 0.;
		}
        template<class _Problem, class _Mesh, class _Result>
        const double fem_solver<_Problem, _Mesh, _Result>::get_value(const _Method2* method2, const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
        {
            if (method2 != nullptr)
                return method2->GetSolution(mesh, res, p);
            return 0.;
        }
        template<class _Problem, class _Mesh, class _Result>
        const double fem_solver<_Problem, _Mesh, _Result>::get_value(const _Method* method2, const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
        {
            if (method2 != nullptr)
                return method2->GetSolution(mesh, res, p);
            return 0.;
        }

        template<class _Problem, class _Mesh, class _Result>
        const double fem_solver<_Problem, _Mesh, _Result>::get_value(const _Mesh& mesh, const _Result& res, const Mesh::Point& p, const int i) const
        {
            if (m_method2 != nullptr)
                return m_method2->GetSolution(mesh, res, p, i);
            return 0.;
        }

		template<class _Problem, class _Mesh, class _Result>
		const Mesh::Point fem_solver<_Problem, _Mesh, _Result>::get_gradvalue(const _Mesh& mesh, const _Result& res, const Mesh::Point& p) const
		{
            if (m_method2 != nullptr)
                return m_method2->GetGradSolution(mesh, res, p);
			return Mesh::Point(0, 0, 0);
		}

        template<class _Problem, class _Mesh, class _Result>
        const Mesh::Point fem_solver<_Problem, _Mesh, _Result>::get_gradvalue(const _Mesh& mesh, const _Result& res, const Mesh::Point& p, const int i) const
        {
            if (m_method2 != nullptr)
                return m_method2->GetGradSolution(mesh, res, p, i);
            return Mesh::Point(0, 0, 0);
        }
	}
}
#endif // !CORENC_SOLVERS_FEM_SOLVER_H_
