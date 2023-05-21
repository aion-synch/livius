#pragma once
#ifndef CORENC_TEST_CASE_REGULAR_MESH_H_
#define CORENC_TEST_CASE_REGULAR_MESH_H_

namespace corenc
{
	namespace tests
	{
		class test_case_regular_mesh
		{
		public:
			test_case_regular_mesh();
			~test_case_regular_mesh();
			const int					construct_mesh() const;
		};
	}
}

#endif // !CORENC_TEST_CASE_REGULAR_MESH_H_
