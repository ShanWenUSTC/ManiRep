#include "ManifoldRepresentation.h"
#include <iostream>
#include <fstream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>
#include <CGAL/Subdivision_method_3.h>

typedef Polyhedron::Halfedge_around_vertex_circulator
	Halfedge_vertex_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator
	Halfedge_face_circulator;
typedef Polyhedron::Vertex_iterator			Vertex_iterator;
typedef Polyhedron::Face_iterator			Face_iterator;
typedef Polyhedron::Halfedge_iterator		Edge_iterator;
typedef Polyhedron::Halfedge_handle			Edge_handle;
typedef Kernel::Point_3						Point;
typedef Kernel::Vector_3					Vector;
typedef Kernel::FT							FT;

using namespace std;
using namespace CGAL;

ManifoldRepresentation::ManifoldRepresentation(void)
{
}


ManifoldRepresentation::~ManifoldRepresentation(void)
{
}

void ManifoldRepresentation::GetModel()
{
	if (charts_.size() != 0)
	{
		charts_.clear();
	}

	charts_.resize(mesh_.size_of_vertices());

	if (manifold_faces_.size() != 0)
	{
		charts_.clear();
	}

	manifold_faces_.resize(mesh_.size_of_facets());

	int index = 0;
	for (Vertex_handle vh = mesh_.vertices_begin(); vh != mesh_.vertices_end();
		++vh)
	{
		charts_[index] = new Chart(vh);
		index++;
	}

	index = 0;
	for (Face_handle fh = mesh_.facets_begin(); fh != mesh_.facets_end();
		++fh)
	{
		manifold_faces_[index] = new ManifoldFace(fh);
		index++;
	}
}

void ManifoldRepresentation::ReadMesh(char* filename)
{
	std::ifstream filein("simple_cube4.off");

	if (!filein)
	{
		cout<<"Cannot open file."<<endl;
	}

	filein >> mesh_; // read the .off

	std::cout<<"read mesh done"<<std::endl;

}

void ManifoldRepresentation::Init()
{

}

void ManifoldRepresentation::CatmullClarkSubdivision()
{
	Polyhedron P(mesh_);
	Subdivision_method_3::CatmullClark_subdivision(P,1);

	Vertex_iterator viter_begin = mesh_.vertices_begin();
	int index=0;
	for (Vertex_iterator viter = mesh_.vertices_begin(); viter != mesh_.vertices_end();
		++viter)
	{
		Point face_point(0, 0, 0);
		Point ptmp(0, 0, 0);
		Halfedge_vertex_circulator hcir = viter->vertex_begin();

		int n = circulator_size(hcir);
		//	cout<<"n "<<n<<endl;

		for (int i=0; i<n; i++, ++hcir)
		{
			Point ep1 = hcir->vertex()->point();
			Point ep2 = hcir->opposite()->vertex()->point();
			point3 edge_point( (ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);
			point2 epoint_local = charts_[index]->CalculateLocalCoordinate(0.5, 0, i);
			//	cout<<"elocal coordinate: "<<epoint_local[0]<<" "<<epoint_local[1]<<endl;
			ChartPoint cpe(epoint_local[0], epoint_local[1], edge_point);
			//	cout<<"edge point: "<<edge_point[0]<<' '<<edge_point[1]<<' '<<edge_point[2]<<endl;
			charts_[index]->points_subdivision_.push_back(cpe);

			Face_handle fh = hcir->facet();
			Halfedge_face_circulator fcir = fh->facet_begin();
			do 
			{ 
				ptmp = ptmp + ( fcir->vertex()->point()-CGAL::ORIGIN);
				// 			face_point[0] += ptmp[0]/4;
				// 			face_point[1] += ptmp[1]/4;
				// 			face_point[2] += ptmp[2]/4;
			} while (++fcir != fh->facet_begin());
			//			face_point = CGAL::ORIGIN + (ptmp - CGAL::ORIGIN)/FT(4);
			face_point = ptmp;
			//	cout<<index<<" face point "<<face_point<<endl;

			point3 epoint(face_point[0]/4, face_point[1]/4, face_point[2]/4);
			point2 point_local = charts_[index]->CalculateLocalCoordinate(0.5, 0.5, i);

			//	cout<<"local coordinate: "<<point_local[0]<<" "<<point_local[1]<<endl;
			ChartPoint cp1(point_local[0], point_local[1], epoint);
			charts_[index]->points_subdivision_.push_back(cp1);
		}

		Vertex_iterator viter2 = P.vertices_begin();
		point3 p_(viter2->point()[0], viter2->point()[1], viter2->point()[2]);
		ChartPoint cptmp(0, 0, p_);
		charts_[index]->points_subdivision_.push_back(cptmp);
		charts_[index]->CalculatePolynomial();
		index++;	
	}
	// 	Face_iterator fiter1 = P.facets_begin();
	// 	Face_iterator fiter2 = fiter1+6;
	// 	int index = fiter2-fiter1;
	// 	cout<<index<<endl;
	//	Subdivision_method_3::CatmullClark_subdivision(P,1);
	// 	for (Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); ++iter)
	// 	{
	// 	}

	// 	for (size_t i=0; i<charts_.size(); i++)
	// 	{
	// 		// face node
	// // 		charts_[i]->vertex_index_->
	// // 		charts_[i]->points_subdivision_.push_back();
	// 
	//	Halfedge_vertex_circulator hcir = charts_[0]->vertex_index_->vertex_begin();
	//	Halfedge_vertex_circulator hcir2 = hcir->opposite();

	//	std::cout<< hcir->vertex()->point()<<std::endl;
	// 		Face_handle fh = hcir->facet();
	// 
	// 		Halfedge_face_circulator fcir = fh->facet_begin();
	// 
	// 		do 
	// 		{
	// 			std::cout<<fcir->vertex()->point()<<std::endl;
	// 		} while (++fcir != fh->facet_begin());
	// 		std::cout<< fh->facet_begin()->vertex()->point()<<std::endl;
	// 		hcir++;
	// 		fh = hcir->facet();
	// 		std::cout<< fh->facet_begin()->vertex()->point()<<std::endl;
	// 		std::cout<< hcir->opposite()->vertex()->point()<<std::endl;

	std::cout<< P;
	// 
	// // 		do 
	// // 				{
	// // 					std::cout<<i<<" "<< hcir->vertex()->point()<<std::endl;
	// // 					hcir2 = hcir->opposite();
	// // 					std::cout<<i<<" "<< hcir2->vertex()->point()<<std::endl;
	// // 				} while (++hcir != charts_[i]->vertex_index_->vertex_begin());
	// // 		
	// // 				Edge_handle eh = charts_[i]->vertex_index_->edge_begin();
	// 				
	// 	}

	// 	Halfedge_vertex_circulator hcir = charts_[0]->vertex_index_->vertex_begin();
	// 	
	// 	std::cout<<hcir->opposite()->vertex()->point()<<std::endl;
	// 	hcir++;
	// 	std::cout<<hcir->opposite()->vertex()->point()<<std::endl;

}

void ManifoldRepresentation::CatmullClarkSubdivision2()
{
	Polyhedron P1(mesh_);
	Polyhedron P2(mesh_);
	Subdivision_method_3::CatmullClark_subdivision(P1,1);
	Subdivision_method_3::CatmullClark_subdivision(P2, 2);

	int num_vertices0 = mesh_.size_of_vertices();
	int num_edge0 = mesh_.size_of_halfedges()/2;
	int num_facets0 = mesh_.size_of_facets();
	int num_vertices1 = P1.size_of_vertices();
	int num_edge1 = P1.size_of_halfedges()/2;
	int num_facets1 = P1.size_of_facets();
	double u,v;

	Vertex_iterator viter = mesh_.vertices_begin();
	for (int i=0; i<charts_.size(); i++)
	{
// 		charts_[i]->subindice_.push_back(i);
// 		charts_[i]->u_.push_back(0);
// 		charts_[i]->v_.push_back(0);

		Halfedge_vertex_circulator vcir = (viter+i)->vertex_begin();
		int n = circulator_size(vcir);

		for (int j=0; j<n; j++, ++vcir)
		{
			// subdivition points on first level
			// face points
			int face_index = vcir->facet()-mesh_.facets_begin();
			Vertex_iterator face_node_iter_inP1 = P1.vertices_begin()+num_vertices0+num_edge0+face_index;
			point3 pface(face_node_iter_inP1->point()[0], face_node_iter_inP1->point()[1], face_node_iter_inP1->point()[2]);
	//		cout<<"face point: "<<pface[0]<<' '<<pface[1]<<' '<<pface[2]<<endl;
			point2 point_local = charts_[i]->CalculateLocalCoordinate(0.5, 0.5, j);
			ChartPoint cpface(point_local[0], point_local[1], pface);
			charts_[j]->points_subdivision_.push_back(cpface);

			// edge points
			int edge_index = (vcir-mesh_.edges_begin())/2;
			Vertex_iterator edge_node_iter_inP1 = P1.vertices_begin()+num_vertices0+edge_index;
			point3 pedge(edge_node_iter_inP1->point()[0], edge_node_iter_inP1->point()[1], edge_node_iter_inP1->point()[2]);
	//		cout<<"edge point: "<<pedge[0]<<' '<<pedge[1]<<' '<<pedge[2]<<endl;
			point_local = charts_[i]->CalculateLocalCoordinate(0, 0.5, j);
			ChartPoint cpedge(point_local[0], point_local[1], pedge);
			charts_[j]->points_subdivision_.push_back(cpedge);

			// subdivision points on second level
			// face points			
 			Halfedge_vertex_circulator vcir_second = face_node_iter_inP1->vertex_begin();
			int n2 = circulator_size(vcir_second);
		//	int testindex = 0;
			while (vcir_second->opposite()->vertex() != edge_node_iter_inP1 )
			{
			//	testindex++;
				++vcir_second;
			}
		//	cout<<"testindex: "<<testindex<<endl;
			int face_index2 = vcir_second->facet()-P1.facets_begin();
			Vertex_iterator face_in_P2 = P2.vertices_begin()+num_vertices1+num_edge1+face_index2;
			point3 pface2(face_in_P2->point()[0], face_in_P2->point()[1], face_in_P2->point()[2]);
			cout<<"face point: "<<pface2[0]<<' '<<pface2[1]<<' '<<pface2[2]<<endl;
			point2 point_loca2 = charts_[i]->CalculateLocalCoordinate(0.25, 0.75, j);
			ChartPoint cpface2(point_loca2[0], point_loca2[1], pface2);
			charts_[j]->points_subdivision_.push_back(cpface2);

			vcir_second++;
			face_index2 = vcir_second->facet()-P1.facets_begin();
			face_in_P2 = P2.vertices_begin()+num_vertices1+num_edge1+face_index2;
			point3 pface2_1(face_in_P2->point()[0], face_in_P2->point()[1], face_in_P2->point()[2]);
			cout<<"face point: "<<pface2_1[0]<<' '<<pface2_1[1]<<' '<<pface2_1[2]<<endl;
			point_loca2 = charts_[i]->CalculateLocalCoordinate(0.25, 0.25, j);
			ChartPoint cpface2_1(point_loca2[0], point_loca2[1], pface2_1);
			charts_[j]->points_subdivision_.push_back(cpface2_1);

			vcir_second++;
			face_index2 = vcir_second->facet()-P1.facets_begin();
			face_in_P2 = P2.vertices_begin()+num_vertices1+num_edge1+face_index2;
			point3 pface2_2(face_in_P2->point()[0], face_in_P2->point()[1], face_in_P2->point()[2]);
			cout<<"face point: "<<pface2_2[0]<<' '<<pface2_2[1]<<' '<<pface2_2[2]<<endl;
			point_loca2 = charts_[i]->CalculateLocalCoordinate(0.75, 0.25, j);
			ChartPoint cpface2_2(point_loca2[0], point_loca2[1], pface2_2);
			charts_[j]->points_subdivision_.push_back(cpface2_2);

			vcir_second++;
			face_index2 = vcir_second->facet()-P1.facets_begin();
			face_in_P2 = P2.vertices_begin()+num_vertices1+num_edge1+face_index2;
			point3 pface2_3(face_in_P2->point()[0], face_in_P2->point()[1], face_in_P2->point()[2]);
			cout<<"face point: "<<pface2_3[0]<<' '<<pface2_3[1]<<' '<<pface2_3[2]<<endl;
			point_loca2 = charts_[i]->CalculateLocalCoordinate(0.75, 0.75, j);
			ChartPoint cpface2_3(point_loca2[0], point_loca2[1], pface2_3);
			charts_[j]->points_subdivision_.push_back(cpface2_3);
			
// 			int n_second = circulator_size(vcir_second);
// 			while (vcir->opposite->vertex() != )
// 			{
// 			}
// 			for (int k=0; k<n; k++, ++vcir_second)
// 			{
// 
// 			}
		}
	}
}

void ManifoldRepresentation::CatmullClarkSubdivision3()
{
	Polyhedron P1(mesh_);
	Polyhedron P2(mesh_);
	Subdivision_method_3::CatmullClark_subdivision(P1, 1);
	Subdivision_method_3::CatmullClark_subdivision(P2, 2);

	int num_vertices0 = mesh_.size_of_vertices();
	int num_edge0 = mesh_.size_of_halfedges()/2;
	int num_facets0 = mesh_.size_of_facets();
	int num_vertices1 = P1.size_of_vertices();
	int num_edge1 = P1.size_of_halfedges()/2;
	int num_facets1 = P1.size_of_facets();
	int num_tmp1 = P2.size_of_vertices();
	int num_tmp2 = P2.size_of_halfedges()/2;
	int num_tmp3 = P2.size_of_facets();
	double u,v;

	Vertex_iterator viter_begin = mesh_.vertices_begin();
	for (int i=0; i<charts_.size(); i++)
	{
		sub_indice subpoint_vertex = {0, 0, i, 0};
		charts_[i]->subindice_.push_back(subpoint_vertex);

		Halfedge_vertex_circulator vcir = (viter_begin+i)->vertex_begin();
		int n = circulator_size(vcir);

		for (int j=0; j<n; j++, ++vcir)
		{
			// subdivision points in first level
			// facet node
			int findex_in_mesh = vcir->facet()-mesh_.facets_begin();
			int findex_in_P1 = num_vertices0+num_edge0+findex_in_mesh;
			sub_indice subpoint_face = {0.5, 0.5, findex_in_P1, 2-j};
			charts_[i]->subindice_.push_back(subpoint_face);

			// test code
// 			Vertex_iterator viter_test = P1.vertices_begin()+findex_in_P1;
// 			cout<< viter_test->point()<<endl;
			
			// edge node
			int eindex_in_mesh = (vcir-mesh_.edges_begin())/2;
			int eindex_in_P1 = num_vertices0+eindex_in_mesh;
			sub_indice subpoint_edge = {0, 0.5, eindex_in_P1, 2-j};
			charts_[i]->subindice_.push_back(subpoint_edge);

			// test code
// 			Vertex_iterator viter_test2 = P1.vertices_begin()+eindex_in_P1;
// 			cout<< viter_test2->point()<<endl;

			// subdivision points on second level
			Halfedge_vertex_circulator vcir_second = (P1.vertices_begin()+findex_in_P1)->vertex_begin();
			int n2 = circulator_size(vcir_second);

			//int testindex = 0;
			while (vcir_second->opposite()->vertex() != P1.vertices_begin()+eindex_in_P1 )
			{
				//testindex++;
				++vcir_second;
			}

			Edge_iterator vtmp = vcir_second->next()->next()->next();
			int eindex_tmp = (vtmp-P1.edges_begin())/2;
			int eindex_tmp_P2 = num_vertices1+eindex_tmp;
			sub_indice subtmp = {0, 0.75, eindex_tmp_P2, 2-j};
			charts_[i]->subindice_.push_back(subtmp);

			vtmp = vcir_second->opposite()->next();
			eindex_tmp = (vtmp-P1.edges_begin())/2;
			eindex_tmp_P2 = num_vertices1+eindex_tmp;
			subtmp.u = 0; subtmp.v = 0.25; subtmp.index_in_P = eindex_tmp_P2;
			charts_[i]->subindice_.push_back(subtmp);

		//	cout<<"test edge: "<<vcir_second->vertex()->point()<<endl;
		//	cout<<"vtmp: "<<vtmp->opposite()->vertex()->point()<<endl;
			// face and edge node on second level
			// first face
			int findex_in_P1_b = vcir_second->facet()-P1.facets_begin();
			int findex_in_P2 = num_vertices1+num_edge1+findex_in_P1_b;
			subpoint_face.u = 0.25; subpoint_face.v = 0.75; subpoint_face.index_in_P = findex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_face);

			// first edge
			int eindex_in_P1_b = (vcir_second-P1.edges_begin())/2;
			int eindex_in_P2 = num_vertices1+eindex_in_P1_b;
			subpoint_edge.u = 0.25; subpoint_edge.v = 0.5; subpoint_edge.index_in_P = eindex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_edge);
			vcir_second++;

			// second face
			findex_in_P1_b = vcir_second->facet()-P1.facets_begin();
			findex_in_P2 = num_vertices1+num_edge1+findex_in_P1_b;
			subpoint_face.u = 0.75; subpoint_face.v = 0.75; subpoint_face.index_in_P = findex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_face);

			// second edge
			eindex_in_P1_b = (vcir_second-P1.edges_begin())/2;
			eindex_in_P2 = num_vertices1+eindex_in_P1_b;
			subpoint_edge.u = 0.5; subpoint_edge.v = 0.75; subpoint_edge.index_in_P = eindex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_edge);
			vcir_second++;

			// third face
			findex_in_P1_b = vcir_second->facet()-P1.facets_begin();
			findex_in_P2 = num_vertices1+num_edge1+findex_in_P1_b;
			subpoint_face.u = 0.75; subpoint_face.v = 0.25; subpoint_face.index_in_P = findex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_face);

			// third edge
			eindex_in_P1_b = (vcir_second-P1.edges_begin())/2;
			eindex_in_P2 = num_vertices1+eindex_in_P1_b;
			subpoint_edge.u = 0.75; subpoint_edge.v = 0.5; subpoint_edge.index_in_P = eindex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_edge);
			vcir_second++;

			// fourth face
			findex_in_P1_b = vcir_second->facet()-P1.facets_begin();
			findex_in_P2 = num_vertices1+num_edge1+findex_in_P1_b;
			subpoint_face.u = 0.25; subpoint_face.v = 0.25; subpoint_face.index_in_P = findex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_face);

			// fourth edge
			eindex_in_P1_b = (vcir_second-P1.edges_begin())/2;
			eindex_in_P2 = num_vertices1+eindex_in_P1_b;
			subpoint_edge.u = 0.5; subpoint_edge.v = 0.25; subpoint_edge.index_in_P = eindex_in_P2;
			charts_[i]->subindice_.push_back(subpoint_edge);
			vcir_second++;

			// subdivision points on origin edge
		//	cout <<testindex <<endl;
		}
	}
}

void ManifoldRepresentation::GetCoordinate(double u, double v)
{
	for (Face_iterator fiter = mesh_.facets_begin(); fiter != mesh_.facets_end(); 
		++fiter)
	{
		Halfedge_face_circulator fcir = fiter->facet_begin();
		point3 coordinate(0, 0, 0);
		double blendsign = 0;
		do 
		{
			Vertex_iterator central_point_iter = fcir->vertex();
			int pindex = central_point_iter-mesh_.vertices_begin();
			//	cout<<"pindex: "<<pindex<<endl;

			// search maniface's indice in charts
			Halfedge_vertex_circulator pcir = central_point_iter->vertex_begin();
			int nsize = circulator_size(pcir);
			int indice;
			Vertex_iterator pitertmp2 = fcir->opposite()->vertex();
			for (int i=0; i<nsize; i++, ++pcir)
			{
				Vertex_iterator pitertmp = pcir->opposite()->vertex();
				if (pitertmp == pitertmp2)
				{
					indice = 2-i;
				//	indice = i;
					
	

					break;
				}
			}
			// 			cout<<"indice "<<indice<<endl;
			// 			cout<<"pindex "<<pindex<<endl;
			// 			cout<<"u v: "<<u<<' '<<v<<endl;
			// 			double ddtmp1 = charts_[pindex]->FuncBlend(u, v);
			// 			cout<<"blend: "<<ddtmp1<<endl;
			// 			point3 ptmp = charts_[pindex]->FuncGeometry(u, v, indice);
			// 			cout<<"geometry: "<<ptmp[0]<<' '<<ptmp[1]<<' '<<ptmp[2]<<endl;
			coordinate += charts_[pindex]->FuncBlend(u, v)* charts_[pindex]->FuncGeometry(u, v, indice);
			blendsign += charts_[pindex]->FuncBlend(u,v);
		//	cout<<"model points: "<<coordinate[0]<<' '<<coordinate[1]<<' '<<coordinate[2]<<endl;
			
			double dtmp = u;
			u = v;
			v = 1-dtmp;

		} while (++fcir != fiter->facet_begin());	
		if(blendsign < 1)
		{
			cout<<"blend wrong"<<endl;
		}
		/*cout<<"blendsign: "<<blendsign<<endl;*/
		model_points_.push_back(coordinate);
	}
}

void ManifoldRepresentation::WriteModel()
{
	double u=0, v=0;
	int n=10;
	for (int i=0; i<n; i++)
	{
		u = 1.0/(n-1)*i;
		for (int j=0; j<n; j++)
		{
			v = 1.0/(n-1)*j;
			GetCoordinate(u, v);
		}
	}

	ofstream fout;
	fout.open("result_test.obj");

	if (!fout)
	{
		cout<<"cannot open file"<<endl;
	}
	for (int i=0; i<model_points_.size(); i++)
	{
		fout<<"v "<<model_points_[i][0]<<" "<<model_points_[i][1]<<" "<<model_points_[i][2]<<endl;
		//	fout<<"v "<<model_points_[i]<<endl;
		//	cout << 'v'<<' '<<model_points_[0]<<' '<<model_points_[1]<<' '<<model_points_[2]<<endl;
	}
	fout.close();

}

void ManifoldRepresentation::test()
{
	// 	for (int i=0; i<charts_.size(); i++)
	// 	{
	// 		for (int j=0; j<charts_[i]->points_subdivision_.size(); j++)
	// 		{
	// 			point3 p = charts_[i]->points_subdivision_[j].point_coordinate;
	// 			cout<<"see: "<<charts_[i]->points_subdivision_[j].x<<' '<<charts_[i]->points_subdivision_[j].y<<' '
	// 				<<p[0]<<' '<<p[1]<<' '<<p[2]<<endl;
	// 		}
	// 		
	// 	}
	// 	for (size_t i=0; i<charts_.size(); i++)
	// 	{
	// 		charts_[i]->CalculatePolynomial();
	// 	}

	// 	ofstream fout("face_point.obj");
	// 	int n = 10;
	// 	double u=0, v=0;
	// 	double d=0;
	// 	for (int i=0; i<n; i++)
	// 	{
	// 		u = 1.0/(n-1)*i;
	// 		for (int j=0; j<n; j++)
	// 		{
	// 			v = 1.0/(n-1)*j;
	// 			point2 p = charts_[0]->CalculateLocalCoordinate(u, v, 0);
	// 			fout<<"v "<<p[0]<<' '<<p[1]<<' '<<d<<endl;
	// 		}
	// 	}

	ofstream fout("single_facet_test.obj");
	int n = 10;
	double u=0, v=0;
	double d=0;
	for (int i=0; i<n; i++)
	{
		u = 1.0/(n-1)*i;
		for (int j=0; j<n; j++)
		{
			v = 1.0/(n-1)*j;

			Face_iterator fiter = mesh_.facets_begin()+2;
			Halfedge_face_circulator fcir = fiter->facet_begin();
			point3 coordinate(0, 0, 0);
			do 
			{
				Vertex_iterator central_point_iter = fcir->vertex();
				int pindex = central_point_iter-mesh_.vertices_begin();
				//	cout<<"pindex: "<<pindex<<endl;

				// search maniface's indice in charts
				Halfedge_vertex_circulator pcir = central_point_iter->vertex_begin();
				int nsize = circulator_size(pcir);
				int indice;
				Vertex_iterator pitertmp2 = fcir->opposite()->vertex();

				for (int k=0; k<nsize; k++, ++pcir)
				{
					Vertex_iterator pitertmp = pcir->opposite()->vertex();
					if (pitertmp == pitertmp2)
					{
						indice = 2-k;
					}
				}
			//	coordinate += charts_[pindex]->FuncBlend(u, v)* charts_[pindex]->FuncGeometry(u, v, indice);
				coordinate += charts_[pindex]->FuncGeometry(u, v, indice);
				
				double dtmp = u;
				u = v;
				v = 1-dtmp;
			}while(++fcir != fiter->facet_begin());
			fout<<"v "<<coordinate[0]<<' '<<coordinate[1]<<' '<<coordinate[2]<<endl;
		}
	}
}

void ManifoldRepresentation::test2()
{
	CatmullClarkSubdivision();

	ofstream fout, foutp;
	fout.open("cc_test.obj");
	foutp.open("cc_test_plane.obj");

	double d1,d2;
	point3 p;
	for (int i=0; i<charts_[0]->points_subdivision_.size(); i++)
	{
		d1 = charts_[0]->points_subdivision_[i].x;
		d2 = charts_[0]->points_subdivision_[i].y;
		p = charts_[0]->points_subdivision_[i].point_coordinate;

		foutp<<"v "<<d1<<' '<<d2<<' '<<endl;
		fout<<"v "<<p[0]<<' '<<p[1]<<' '<<p[2]<<endl;
	}
	fout.close();
	foutp.close();
}

void ManifoldRepresentation::test3()
{
// 	CatmullClarkSubdivision();
 	ofstream fout;
 	fout.open("single_face_in_single_chart.obj");

	int n=10;
	double u,v;
	for (int i=0; i<n; i++)
	{
		u = 1.0/(n-1)*i;
		for (int j=0; j<n; j++)
		{
			v = 1.0/(n-1)*j;
			point3 p = charts_[0]->FuncGeometry(u, v, 0);
			fout<<"v "<<p[0]<<' '<<p[1]<<' '<<p[2]<<endl;
		}
	}
	fout.close();
// 	Polyhedron P(mesh_);
// 	Subdivision_method_3::CatmullClark_subdivision(P,1);
// 
// 	cout<<"mesh: "<<mesh_<<endl;
// 	cout<<"P: "<<P<<endl;
}

void ManifoldRepresentation::test4()
{
	Polyhedron P(mesh_);
	Subdivision_method_3::CatmullClark_subdivision(P, 2);
	ofstream fout("sub_test.obj");
	ofstream fout2("sub_test_3d.obj");

	double dtmp=0;
	double u, v;
	int findex, meshindex;
	int n = charts_.size();
	
	for (int i=0; i<n; i++)
	{
		charts_[i]->points_subdivision_.clear();
		for (int j=0; j<charts_[i]->subindice_.size(); j++)
		{
			u = charts_[i]->subindice_[j].u;
			v = charts_[i]->subindice_[j].v;
			findex = charts_[i]->subindice_[j].findex;
			meshindex = charts_[i]->subindice_[j].index_in_P;
		//	cout<<"meshindex: "<<meshindex<<endl;

			point2 p = charts_[i]->CalculateLocalCoordinate(u, v, findex);
	//		fout<<"v "<<p[0]<<' '<<p[1]<<' '<<dtmp<<endl;

			Point CGp = (P.vertices_begin()+meshindex)->point();
			point3 p3(CGp.x(), CGp.y(), CGp.z());
	//		fout2<<"v "<<p3[0]<<' '<<p3[1]<<' '<<p3[2]<<endl;
			ChartPoint cptmp(p[0], p[1], p3);
			charts_[i]->points_subdivision_.push_back(cptmp);
		}

		charts_[i]->CalculatePolynomial();
	}

	fout.close();
	fout2.close();
}

void ManifoldRepresentation::test5()
{
	ofstream fout("single chart approximate.obj");

	int index = 7;
	double u=0,v =0;
	int n=10;
	for (int i=0; i<n; i++)
	{
		u = 1.0*i/(n-1);
		for (int j=0; j<n; j++)
		{
			v = 1.0*j/(n-1);
			for (int k=0; k<3; k++)
			{
				point3 p3 = charts_[index]->FuncGeometry(u, v, k);
				fout<<"v "<<p3[0]<<' '<<p3[1]<<' '<<p3[2]<<endl;
			}
		}
	}
	
	fout.close();

	ofstream fout2("cc_approximate.obj");
	
	int vindex=1; 
	for(int k=0; k<charts_[vindex]->points_subdivision_.size(); k++)
	{
		ChartPoint cptmp = charts_[vindex]->points_subdivision_[k];
		double dx = cptmp.x;
		double dy = cptmp.y;
		
		int index = 0;
		double x=0, y=0, z=0;
		for (int i=0; i<4+1; i++)
		{
			for (int j=0; j<4+1-i; j++)
			{
				x += charts_[vindex]->basis_coefficients_x_[index]*pow(dx, i)*pow(dy, j);
				y += charts_[vindex]->basis_coefficients_y_[index]*pow(dx, i)*pow(dy, j);
				z += charts_[vindex]->basis_coefficients_z_[index]*pow(dx, i)*pow(dy, j);
				index++;
			}
		}

		fout2<<"v "<<x<<' '<<y<<' '<<z<<endl;
		
	}
	fout2.close();
}

void ManifoldRepresentation::test6()
{
	cout<<"mesh: "<<mesh_<<endl;

	Vertex_iterator viter = mesh_.vertices_begin()+5;

	Halfedge_vertex_circulator vcir = viter->vertex_begin();
	int n = circulator_size(vcir);

	for (int i=0; i<n; i++, ++vcir)
	{
		cout<<"test6: "<<vcir->opposite()->vertex()->point()<<endl;
	}
}

void ManifoldRepresentation::test7()
{
	int n=10;
	double v = 0, u;
	for (int i=0; i<n-1; i++)
	{
		u = 1.0*i/(n-1);
		point2 p = charts_[0]->CalculateLocalCoordinate(u, v, 0);
		cout<<"test7: "<<u<<' '<<v<<' '<<p<<endl;
	}

// 	int n=30;
// 	double u, v;
// 	ofstream fout("test7.obj");
// 	for (int i=0; i<n; i++)
// 	{
// 		u = i*1.0/(n-1);
// 		for (int j=0; j<n; j++)
// 		{
// 			v = j*1.0/(n-1);
// 			Face_iterator fiter = mesh_.facets_begin()+5;
// 
// 			Halfedge_face_circulator fcir = fiter->facet_begin();
// 			point3 coordinate(0, 0, 0);
// 
// 			do 
// 			{
// 				Vertex_iterator central_point_iter = fcir->vertex();
// 				int pindex = central_point_iter-mesh_.vertices_begin();
// 				//	cout<<"pindex: "<<pindex<<endl;
// 
// 				// search maniface's indice in charts
// 				Halfedge_vertex_circulator pcir = central_point_iter->vertex_begin();
// 				int nsize = circulator_size(pcir);
// 				int indice;
// 				Vertex_iterator pitertmp2 = fcir->opposite()->vertex();
// 				for (int k=0; k<nsize; k++, ++pcir)
// 				{
// 					Vertex_iterator pitertmp = pcir->opposite()->vertex();
// 					if (pitertmp == pitertmp2)
// 					{
// 						indice = 2-k;
// 					}
// 				}
// 				// 			cout<<"indice "<<indice<<endl;
// 				// 			cout<<"pindex "<<pindex<<endl;
// 				// 			cout<<"u v: "<<u<<' '<<v<<endl;
// 				// 			double ddtmp1 = charts_[pindex]->FuncBlend(u, v);
// 				// 			cout<<"blend: "<<ddtmp1<<endl;
// 				// 			point3 ptmp = charts_[pindex]->FuncGeometry(u, v, indice);
// 				// 			cout<<"geometry: "<<ptmp[0]<<' '<<ptmp[1]<<' '<<ptmp[2]<<endl;
// 				coordinate += charts_[pindex]->FuncBlend(u, v)* charts_[pindex]->FuncGeometry(u, v, indice);
// 				//	cout<<"model points: "<<coordinate[0]<<' '<<coordinate[1]<<' '<<coordinate[2]<<endl;
// 
// 				double dtmp = u;
// 				u = v;
// 				v = 1-dtmp;
// 
// 			} while (++fcir != fiter->facet_begin());	
// 			fout<<"v "<<coordinate[0]<<' '<<coordinate[1]<<' '<<coordinate[2]<<endl;
// 		}
// 	}	
// 
// 	fout.close();
}

void ManifoldRepresentation::test8()
{
	Polyhedron P1(mesh_);
	Polyhedron P2(mesh_);
	Subdivision_method_3::CatmullClark_subdivision(P1, 1);
	Subdivision_method_3::CatmullClark_subdivision(P2, 2);
	cout<<"mesh: "<<mesh_<<endl;
	cout<<"P1: "<<P1<<endl;
	cout<<"P2: "<<P2<<endl;
// 	ofstream fout("sub_test.obj");
// 	ofstream fout2("sub_test_3d.obj");
// 
// 	double dtmp=0;
// 	double u, v;
// 	int findex, meshindex;
// 	int n = charts_.size();
// 
// 	for (int i=2; i<3; i++)
// 	{
// 		charts_[i]->points_subdivision_.clear();
// 		for (int j=0; j<charts_[i]->subindice_.size(); j++)
// 		{
// 			u = charts_[i]->subindice_[j].u;
// 			v = charts_[i]->subindice_[j].v;
// 			findex = charts_[i]->subindice_[j].findex;
// 			meshindex = charts_[i]->subindice_[j].index_in_P;
// 			//	cout<<"meshindex: "<<meshindex<<endl;
// 
// 			point2 p = charts_[i]->CalculateLocalCoordinate(u, v, findex);
// 			//		fout<<"v "<<p[0]<<' '<<p[1]<<' '<<dtmp<<endl;
// 
// 			Point CGp = (P.vertices_begin()+meshindex)->point();
// 			point3 p3(CGp.x(), CGp.y(), CGp.z());
// 					fout2<<"v "<<p3[0]<<' '<<p3[1]<<' '<<p3[2]<<endl;
// 			ChartPoint cptmp(p[0], p[1], p3);
// 			charts_[i]->points_subdivision_.push_back(cptmp);
// 		}
// 
// 		charts_[i]->CalculatePolynomial();
// 	}
// 
// 	fout.close();
// 	fout2.close();
}