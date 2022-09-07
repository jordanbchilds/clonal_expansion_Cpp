#include <cmath>
#include <cassert>
#include <poplar/Vertex.hpp>
#include <algorithm>
#include <array>
#include <ipu_vector_math>
#include <ipudef.h>
#include <ipu_builtins.h>

using namespace poplar;
using namespace std;

class example_vertex : public poplar::Vertex
{
public:
	Input<Vector<float>> x ;
	
    Output<Vector<float>> out ;

	bool compute()
	{
		for(int i=0; i<x.size(); ++i){
			out[i] = theta[i];
		}

		return true;
	}
};
