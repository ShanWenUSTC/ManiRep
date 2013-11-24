// ManiRep.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>
#include <iostream>
#include <fstream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
using namespace std;
using namespace CGAL;

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc != 2) 
	{
		cout << "Usage: CatmullClark_subdivision d < filename" << endl;
		cout << " d: the depth of the subdivision (0 < d < 10)" << endl;
		cout << " filename: the input mesh (.off)" << endl;
		return 0;
	}
	int d = argv[1][0] - '0';
	Polyhedron P;
	cin >> P; // read the .off
	//	Subdivision_method_3::CatmullClark_subdivision(P,d);
	//	cout << P; // write the .off

	//	CGAL::ostream out( "test.off");
	std::ofstream file("test.off");
	file << P;

	file.close();
	return 0;
}

