#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Root_of_traits.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::Locate_type Locate_type;
typedef Delaunay::Cell_handle Cell_handle;

typedef K::FT FT;
typedef CGAL::Point_3<K> Point_3;
typedef CGAL::Sphere_3<K> Sphere_3;

bool comp_sph(const Sphere_3 &a, const Sphere_3 &b) {
  return a.squared_radius()>b.squared_radius();
}

int main(int argc, char *argv[])
{
  std::string input=argv[1];
  std::string output=argv[2];
  FT lowRange=0;	// Mpc/h
  FT highRange=2592;	// Mpc/h
  FT cpyRange=100;	// Mpc/h
//  FT outRange=0;
  FT boxsize=highRange-lowRange;
  FT radius;
  FT Vol;

  std::vector<Point> P;
  std::vector<Point> Pin;
  std::vector<Sphere_3> SP;
  Point_3 PP[4];
  Sphere_3 Sph;
  std::stringstream str;
  std::size_t size, i;

  if(argc != 3) {
    std::cout<<"Usage: "<<argv[0]<<" input_catalog output_catalog"<<std::endl;
    return -1;
  }

  std::cout<<"Reading file: "<<input<<std::endl;
  std::ifstream iFile(argv[1],std::ios::in);
  assert(iFile);

  while(iFile>>PP[0])
    P.push_back(PP[0]);
  size=P.size();
  std::cout<<"Number of records: "<<size<<std::endl<<std::endl;
  iFile.close();

  std::cout<<"Duplicating boundaries for periodic condition"<<std::endl;
  for(i=0;i<size;i++)
    if(P[i].x()<lowRange+cpyRange)
      P.push_back(Point_3(P[i].x()+boxsize,P[i].y(),P[i].z()));
  size=P.size();
  for(i=0;i<size;i++) 
    if(P[i].x()>=highRange-cpyRange && P[i].x()<highRange)
      P.push_back(Point_3(P[i].x()-boxsize,P[i].y(),P[i].z()));
  size=P.size();
  for(i=0;i<size;i++) 
    if(P[i].y()<lowRange+cpyRange)
      P.push_back(Point_3(P[i].x(),P[i].y()+boxsize,P[i].z()));
  size=P.size();
  for(i=0;i<size;i++) 
    if(P[i].y()>=highRange-cpyRange && P[i].y()<highRange)
      P.push_back(Point_3(P[i].x(),P[i].y()-boxsize,P[i].z()));
  size=P.size();
  for(i=0;i<size;i++) 
    if(P[i].z()<lowRange+cpyRange)
      P.push_back(Point_3(P[i].x(),P[i].y(),P[i].z()+boxsize));
  size=P.size();
  for(i=0;i<size;i++) 
    if(P[i].z()>=highRange-cpyRange && P[i].z()<highRange)
      P.push_back(Point_3(P[i].x(),P[i].y(),P[i].z()-boxsize));
  size=P.size();
  std::cout<<"Number of records: "<<size<<std::endl<<std::endl;

  std::cout<<"Building Delaunay Triangulation."<<std::endl;
  Delaunay T(P.begin(), P.end());
  assert(T.is_valid());
  P.clear();

  std::cout<<"Number of vertices: "<<T.number_of_vertices()<<std::endl;
  std::cout<<"Number of cells: "<<T.number_of_finite_cells()<<std::endl<<std::endl;;
  std::cout<<"Writing to file "<<output<<std::endl;

  std::ofstream oFile(argv[2],std::ios::out);
  assert(oFile);
  oFile.precision(12);

  for(Delaunay::Finite_cells_iterator cell=T.finite_cells_begin();cell!=T.finite_cells_end();cell++) {
    for(i=0;i<4;i++)
      PP[i] = cell->vertex(i)->point();
    Sph = Sphere_3(PP[0],PP[1],PP[2],PP[3]);
    PP[0]=Sph.center();
    if(PP[0].x()<lowRange || PP[0].x()>=highRange ||
        PP[0].y()<lowRange || PP[0].y()>=highRange ||
        PP[0].z()<lowRange || PP[0].z()>=highRange)
      continue;

//    str<<Sph.center()<<" "<<CGAL::make_sqrt(Sph.squared_radius())<<std::endl;
    oFile<<PP[0]<<" "<<CGAL::sqrt(CGAL::to_double(Sph.squared_radius()))<<std::endl;
//    str<<PP[0]<<" "<<CGAL::make_sqrt(Sph.squared_radius())<<std::endl;
/*    size++;
    if(size==1000000) {
      oFile<<str.str();
      size=0;
      str.str(std::string());
      str.clear();
    }
*/
  }
//  if(size!=0) oFile<<str.str();

  oFile.close();

  std::cout<<"Finished successfully!"<<std::endl;
  return 0;
}

