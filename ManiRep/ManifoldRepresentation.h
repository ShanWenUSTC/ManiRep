#pragma once
#include <vector>
#include "Chart.h"
#include "ManifoldFace.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

typedef CGAL::Simple_cartesian<double>		Kernel;
typedef CGAL::Polyhedron_3< Kernel,	
			CGAL::Polyhedron_items_3,
			CGAL::HalfedgeDS_vector> Polyhedron;

using std::vector;

class ManifoldRepresentation
{
public:
	ManifoldRepresentation(void);
	~ManifoldRepresentation(void);

	void GetModel();
	void ReadMesh(char* filename);
	void Init();
	void CatmullClarkSubdivision();
	void CatmullClarkSubdivision2();
	void CatmullClarkSubdivision3();
	void test();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void GetCoordinate(double u, double v);
	void WriteModel();

private:
	vector<ManifoldFace*> manifold_faces_;
	vector<Chart*> charts_;
	vector<point3> model_points_;
	Polyhedron mesh_;
};

