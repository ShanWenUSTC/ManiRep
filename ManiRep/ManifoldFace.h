#pragma once
#include <vector>
#include <Eigen/Dense>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Simple_cartesian<double>		Kernel;
typedef CGAL::Polyhedron_3< Kernel,	
	CGAL::Polyhedron_items_3,
	CGAL::HalfedgeDS_vector> Polyhedron;
typedef Polyhedron::Face_handle			Face_handle;
typedef Polyhedron::Halfedge_around_facet_circulator
	Halfedge_face_circulator;
using std::vector;
using Eigen::Vector3d;

typedef Vector3d point3;

class Chart;
class ManifoldFace
{
public:
	ManifoldFace(Face_handle fhandle);
	~ManifoldFace(void);

	point3 GetCoordinate(double u, double v);
	void CCsubdivision();

private:
	vector<Chart*> neighbor_charts_;
	vector<int> index_maniface;
	vector<point3> neighbor_points_;
	vector<point3> points_subdivision_;
	Face_handle face_index_;
};

