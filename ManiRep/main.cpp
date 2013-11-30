// ManiRep.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "ManifoldRepresentation.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>
#include <iostream>
#include <fstream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

typedef CGAL::Simple_cartesian<double> Kernel;
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
using namespace std;
using namespace CGAL;
#include <Eigen/Dense>
#include <math.h>
using Eigen::MatrixXd;
using namespace::std;

int _tmain(int argc, _TCHAR* argv[])
{
// 	if (argc != 2) 
// 	{
// 		cout << "Usage: CatmullClark_subdivision d < filename" << endl;
// 		cout << " d: the depth of the subdivision (0 < d < 10)" << endl;
// 		cout << " filename: the input mesh (.off)" << endl;
// 		return 0;
// 	}/
// 	int d = argv[1][0] - '0';
//  	ifstream filein;
//  	filein.open("simple_cube4.off");
// 	if (filein)
// 	{
// 		cout<<"filein OK"<<endl;
// 	}
// 
	/*Polyhedron P;*/

//	cout<<"OK"<<endl;
	/*filein >> P; */// read the .off
//	filein.close();
	/*Subdivision_method_3::CatmullClark_subdivision(P,2);*/
//	cout << P; // write the .off

	//	CGAL::ostream out( "test.off");
// 	std::ofstream file("test.off");
// 	file << P;
// 
// 	file.close();

	ManifoldRepresentation mani_rep;
 	mani_rep.ReadMesh("simple_cube4.obj");
 	mani_rep.GetModel();
 	mani_rep.CatmullClarkSubdivision3();
	mani_rep.test4();
	mani_rep.WriteModel();
//	mani_rep.test();

//	mani_rep.CatmullClarkSubdivision();
//	mani_rep.test2();
//	mani_rep.test3();

//	mani_rep.test3();
//	mani_rep.test5();
//	mani_rep.test6();
//	mani_rep.test8();
//	mani_rep.test7();
// 	MatrixXd m(3,3);
// 	m<<1,2,3,
// 		4,5,6,
// 		7,8,9;
// 	m(2,2) = pow(1.0, 8)*pow(0.777, 9);
// 	cout<<m(2,2)<<endl;
	return 0;
}

