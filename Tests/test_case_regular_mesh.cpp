#include "test_case_regular_mesh.h"
#include "../CoreNCFEM/Grids/RegularMesh.h"
using namespace std;
using namespace corenc;
using namespace tests;
using namespace Mesh;
test_case_regular_mesh::test_case_regular_mesh()
{

}

test_case_regular_mesh::~test_case_regular_mesh()
{

}

const int test_case_regular_mesh::construct_mesh() const
{
	CRegularMesh mesh{ 0, 0, 1, 1, 3, 3 };
	return 0;
}