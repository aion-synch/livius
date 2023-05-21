#ifndef CORENC_FINITESOLVER_H
#define CORENC_FINITESOLVER_H

namespace corenc
{
	template<class Method, class Mesh, class Solver>
	class CFiniteSolver
	{
	public:
		CFiniteSolver() {};
		~CFiniteSolver() {};
		void						Solve();
	private:
		Method*						m_method;
		Mesh*						m_mesh;
		Solver*						m_solver;
	};

	template<class Method, class Mesh, class Solver>
	void CFiniteSolver<Method, Mesh, Solver>::Solve()
	{
		m_method->Assemble();
		return;
	}
}

#endif // !CORENC_FINITESOLVER_H

