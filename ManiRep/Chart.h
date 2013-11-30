#pragma once
#include <vector>
#include <Eigen/Dense>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_3.h>

using std::vector;
using Eigen::VectorXd;

typedef Eigen::Vector3d point3;
typedef Eigen::Vector2d point2;
typedef CGAL::Simple_cartesian<double>		Kernel;
typedef CGAL::Polyhedron_3< Kernel,	
	CGAL::Polyhedron_items_3,
	CGAL::HalfedgeDS_vector> Polyhedron;
typedef Polyhedron::Vertex_handle			Vertex_handle;

class ChartPoint  
{
public:
	ChartPoint(double dx, double dy, point3 p);
	~ChartPoint();
	double x;
	double y;
	point3 point_coordinate;
};

struct sub_indice	 
{
	double u;
	double v;
	int index_in_P;
	int findex;
};

class ManifoldFace;
class Chart
{
public:
	Chart(Vertex_handle vhandle);
	~Chart(void);

	double FuncBlend(double u, double v);
	point3 FuncGeometry(double u, double v, int faceIndex);
	point2 CalculateLocalCoordinate(double u, double v, int faceIndex);
	void CalculatePolynomial();

	Vertex_handle vertex_index_;
	vector<ChartPoint> points_subdivision_;

public:
	double h(double s);

	int valence_;
	int basis_degree_;
	VectorXd basis_coefficients_x_;
	VectorXd basis_coefficients_y_;
	VectorXd basis_coefficients_z_;
	point3 point_central_;
	vector<point3> points_control_;
	vector<ManifoldFace*> faces_;
	vector<sub_indice> subindice_;	
};

