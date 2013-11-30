#include "ManifoldFace.h"
#include <iostream>
#include "Chart.h"

using std::cout;
using std::endl;
using Eigen::Matrix3d;

ManifoldFace::ManifoldFace(Face_handle fhandle)
{
	face_index_ = fhandle;
}

ManifoldFace::~ManifoldFace(void)
{
}

point3 ManifoldFace::GetCoordinate(double u, double v)
{
	if ((u>1)||(v>1))
	{
		cout<<"Wrong in bilinear coordinate!"<<endl;
	}

	point3 p(0,0,0);
	
	Halfedge_face_circulator fcir = face_index_->facet_begin();
	for (size_t i=0; i<neighbor_charts_.size(); i++)
	{
		p += neighbor_charts_[i]->FuncBlend(u, v) 
			* neighbor_charts_[i]->FuncGeometry(u, v, index_maniface[i]);
	}

	return p;
}

void ManifoldFace::CCsubdivision()
{
	Matrix3d iner_point;
}
