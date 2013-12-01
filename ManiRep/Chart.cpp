#include "Chart.h"
#include <math.h>
#include <Eigen/Dense>
#include "ManifoldFace.h"
#include <iostream>

using Eigen::Matrix2d;
using Eigen::MatrixXd;
using namespace std;

#define PI 3.1415926
typedef Matrix2d tranform;

Chart::Chart(Vertex_handle vhandle)
{
	vertex_index_ = vhandle;
	valence_ = 3;
	basis_degree_ = 2;
}


Chart::~Chart(void)
{
}

double Chart::FuncBlend(double u, double v)
{
	double u_ita, v_ita;
	if (u<0.125)
	{
		u_ita = 1;
	}
	else if (u>0.875)
	{
		u_ita = 0;
	}
	else
	{
		double dtmp1 = (u-0.125)/0.75;
	//	double dtmp2 = (1-u-0.125)/0.75;
		double dtmp2 = 1-dtmp1;
		u_ita = h(dtmp1)/(h(dtmp1) + h(dtmp2));
	}

	if (v<0.125)
	{
		v_ita = 1;
	}
	else if (v>0.875)
	{
		v_ita = 0;
	}
	else
	{
		double dtmp1 = (v-0.125)/0.75;
	//	double dtmp2 = (1-v-0.125)/0.75;
		double dtmp2 = 1-dtmp1;
		v_ita = h(dtmp1)/(h(dtmp1) + h(dtmp2));
	}
	
	return u_ita*v_ita;
}

point3 Chart::FuncGeometry(double u, double v, int faceIndex)
{
	point2 p = CalculateLocalCoordinate(u, v, faceIndex);
//	cout<<"ptest: "<<p[0]<<' '<<p[1]<<endl;
	int index = 0;
	double x=0, y=0, z=0;
	for (int i=0; i<basis_degree_+1; i++)
	{
		for (int j=0; j<basis_degree_+1-i; j++)
		{
			x += basis_coefficients_x_[index]*pow(p[0], i)*pow(p[1], j);
			y += basis_coefficients_y_[index]*pow(p[0], i)*pow(p[1], j);
			z += basis_coefficients_z_[index]*pow(p[0], i)*pow(p[1], j);
			index++;
		}
	}

	point3 ptmp(x, y, z);
	return ptmp;
}

point2 Chart::CalculateLocalCoordinate(double u, double v, int faceIndex)
{
	if ((u<0.00000001) && (v<0.00000001))
	{
		point2 p(0,0);
		return p;
	}
	point2 p0, p1, p2, p3;
	p0[0] = 0; p0[1] = 0;
	p1[0] = 1; p1[1] = 0;
	p2[0] = 1; p2[1] = 1;
	p3[0] = 0; p3[1] = 1;
// 	tranform rtmp1, rtmp2;
// 	double ctmp = cos(2*PI/valence_);
// 	double stmp = sin(2*PI/valence_);
// 	double ctmp_half = cos(PI/valence_);
// 	double stmp_half = sin(PI/valence_);
// 
// 	rtmp1 << ctmp, -stmp,
// 		  stmp, ctmp;
// 	rtmp2 << ctmp_half, -stmp_half,
// 		stmp_half, ctmp_half;
// 	p3 = rtmp1*p1;
// 	p2 = 2*ctmp_half*rtmp2*p1;
	/*cout<<"p: "<<p2[0]<<' '<<p2[1]<<' '<<p3[0]<<' '<<p3[1]<<endl;*/
	point2 ptmp = (1-u)*(1-v)*p0 + u*(1-v)*p1 + u*v*p2 + (1-u)*v*p3;
	double dxtmp = ptmp[0];
	double dytmp = ptmp[1];
//	cout<<ptmp<<endl;
//	cout<<"point: "<<ptmp<<endl;

//	tranform rotate, linear;

// 	ctmp = cos(2.0*faceIndex*PI/valence_);
// 	stmp = sin(2.0*faceIndex*PI/valence_);
// 	rotate << ctmp, -stmp,
// 		stmp, ctmp;
//  	linear << cos(PI/4)/cos(PI/valence_), 0,
//  		0, sin(PI/4)/sin(PI/valence_);
	
//	cout<<"linear: "<<linear<<endl;
// 	point2 central = (p0+p1+p2+p3)/4;
// 	point2 ptmp = central+linear*(p_init-central);


// 
// /*	cout<<"p linear: "<<ptmp[0]<<' '<<ptmp[1]<<endl;*/
// 
	double angle;
	double length = sqrt(ptmp[0]*ptmp[0]+ptmp[1]*ptmp[1]);
	if (ptmp[1]>-0.0000001)
	{
		angle = acos(ptmp[0]/length);
	} 
	else
	{
		angle = 2*PI - acos(ptmp[0]/length);
	//	cout<<"here"<<endl;
	}

	double ctmp = cos(angle*4/valence_);
	double stmp = sin(angle*4/valence_);
	length = pow(length, 4/valence_);

	point2 p(ctmp*length, stmp*length);
	ctmp = cos(2.0*faceIndex*PI/valence_);
	stmp = sin(2.0*faceIndex*PI/valence_);
	tranform rotate;
	rotate << ctmp, -stmp,
		stmp, ctmp;
	p = rotate*p;

	return p;
}

void Chart::CalculatePolynomial()
{
// 	vector<point3> points_fitting;
// 
// 	point2 p0, p1, p2, p3;
// 	p0 = CalculateLocalCoordinate(0, 0, 0);
// 	p1 = CalculateLocalCoordinate(1, 0, 0);
// 	p2 = CalculateLocalCoordinate(1, 1, 0);
// 	p3 = CalculateLocalCoordinate(0, 1, 0);
	
	//int num_row = 12*valence_+1;
	int num_row = points_subdivision_.size();
	int num_column = (basis_degree_+1)*(basis_degree_+2)/2;
	MatrixXd A(num_row, num_column);
	VectorXd b_x(num_row);
	VectorXd b_y(num_row);
	VectorXd b_z(num_row);
	basis_coefficients_x_.resize(num_column);
	basis_coefficients_y_.resize(num_column);
	basis_coefficients_z_.resize(num_column);

	for (int i=0; i<num_row; i++)
	{
		int index =0;
	//	cout<<"points_sub x: "<<points_subdivision_[i].x<<endl;
	//	cout<<"points_sub y: "<<points_subdivision_[i].y<<endl;
		for (int j=0; j<basis_degree_+1; j++)
		{
			for (int k=0; k<basis_degree_+1-j; k++)
			{
				A(i, index) = pow(points_subdivision_[i].x, j)
						*pow(points_subdivision_[i].y, k);
				
 			//	cout<<"powx: "<<pow(points_subdivision_[i].x, j)<<endl;
 			//	cout<<"powy: "<<pow(points_subdivision_[i].y, k)<<endl;
				/*cout<<"x*y: "<<pow(points_subdivision_[i].x, j)*pow(points_subdivision_[i].y, k)<<endl;*/
				
	//			cout<<"A: "<<i<<' '<<index<<' '<<A(i, index)<<endl;
				index++; 
			}
		}

		b_x[i] = points_subdivision_[i].point_coordinate[0];
		b_y[i] = points_subdivision_[i].point_coordinate[1];
		b_z[i] = points_subdivision_[i].point_coordinate[2];

// 		b_x[i] = points_subdivision_[i].x;
// 		b_y[i] = points_subdivision_[i].y;
// 		b_z[i] = pow(points_subdivision_[i].x, 4)+pow(points_subdivision_[i].y, 4);

	//	cout<<"b: "<<b_x[i]<<' '<<b_y[i]<<' '<<b_z[i]<<endl;
	}

//  	for (int i=0; i<b_x.size(); i++)
//  	{
//  		cout<<"b: "<<b_x[i]<<' '<<b_y[i]<<' '<<b_z[i]<<endl;
//  	}
	basis_coefficients_x_ = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_x);
	basis_coefficients_y_ = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_y);
	basis_coefficients_z_ = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_z);

	for(size_t i=0; i<basis_coefficients_x_.size(); i++)
	{
		cout<<"bcs: "<<basis_coefficients_x_[i]<<' '<<basis_coefficients_y_[i]<<' '<<basis_coefficients_z_[i]<<endl;
	}
}

double Chart::h(double s)
{
	return exp( 2*exp(-1/s)/(s-1) );
}

ChartPoint::ChartPoint(double dx, double dy, point3 p)
{
	x = dx;
	y = dy;
	point_coordinate = p;
}

ChartPoint::~ChartPoint()
{

}
