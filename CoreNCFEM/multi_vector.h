#pragma once
#ifndef CORENC_MULTI_VECTOR_H_
#define CORENC_MULTI_VECTOR_H_
#include <vector>
#include <cstdarg>
#include <cstddef>
namespace corenc
{
	template<class T>
	class multi_vector
	{
	public:
		multi_vector();
		// dim = 1 vector, dim = 2 matrix, etc
		// block x ... x block; dim times
		multi_vector(const size_t block, const size_t dim);
		multi_vector(const size_t dim);
		~multi_vector();
		const T				get(const size_t i...) const;
		const T				get(const std::vector<size_t>& i) const;
		const int			set(const T& element, const std::vector<size_t>& index);
		//const int			set(const T& element, const size_t i...);
		const int			fill_inc();
		void				resize(const size_t block);
		void				resize(const size_t block, const size_t dim);
		const size_t		size() const;
		const size_t		totalsize() const;
	private:
		std::vector<T>		m_vector;
		size_t				m_dim;
		size_t				m_block;
		size_t				m_totalsize;
	};

	template<class T>
	multi_vector<T>::multi_vector()
	{

	}
	template<class T>
	multi_vector<T>::multi_vector(const size_t block, const size_t dim)
	{
		m_block = block;
		m_dim = dim;
		m_totalsize = 1;
		for (size_t i = 0; i < m_dim; ++i, m_totalsize *= block);
		m_vector.resize(m_totalsize);
	}
	template<class T>
	multi_vector<T>::multi_vector(const size_t dim)
	{
		m_block = 0;
		m_dim = dim;
		m_totalsize = 0;
	}
	template<class T>
	multi_vector<T>::~multi_vector()
	{

	}
	template<class T>
	const size_t multi_vector<T>::size() const
	{
		return m_block;
	}
	template<class T>
	const size_t multi_vector<T>::totalsize() const
	{
		return m_totalsize;
	}
	template<class T>
	const T multi_vector<T>::get(const size_t i...) const
	{
		va_list args;
		va_start(args, i);

		va_end(args);
		return m_vector[i];
	}
	template<class T>
	const T multi_vector<T>::get(const std::vector<size_t>& i) const
	{
		if (i.size() != m_dim)
			return T(0);
		size_t ind = 0;
		for (size_t j = 0; j < m_dim; ++j)
		{
			size_t l = 1;
			const int lim = m_dim - j - 1;
			for (int k = 0; k < lim; ++k, l *= m_block);
			ind += i[j] * l;
		}
		return m_vector[ind];
	}
	template<class T>
	void multi_vector<T>::resize(const size_t block)
	{
		m_block = block;
		m_totalsize = 1;
		for (size_t i = 0; i < m_dim; ++i, m_totalsize *= block);
		m_vector.resize(m_totalsize);
	}
	template<class T>
	void multi_vector<T>::resize(const size_t block, const size_t dim)
	{
		m_block = block;
		m_dim = dim;
		m_totalsize = 1;
		for (size_t i = 0; i < m_dim; ++i, m_totalsize *= block);
		m_vector.resize(m_totalsize);
	}
	template<class T>
	const int multi_vector<T>::fill_inc()
	{
		for (size_t i = 0; i < m_totalsize; ++i)
			m_vector[i] = i;
		return 0;
	}
	template<class T>
	const int multi_vector<T>::set(const T& element, const std::vector<size_t>& i)
	{
		if (i.size() != m_dim)
			return 1;
		size_t ind = 0;
		for (size_t j = 0; j < m_dim; ++j)
		{
			size_t l = 1;
			const int lim = m_dim - j - 1;
			for (int k = 0; k < lim; ++k, l *= m_block);
			ind += i[j] * l;
		}
		m_vector[ind] = element;
		return 0;
	}
}
#endif // !CORENC_MULTI_VECTOR_H_
