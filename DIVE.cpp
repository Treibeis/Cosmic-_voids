#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
//#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/number_utils.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::Locate_type Locate_type;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;

typedef K::FT FT;
typedef CGAL::Point_3<K> Point_3;
typedef CGAL::Sphere_3<K> Sphere_3;
typedef CGAL::Search_traits_3<K> ST;
typedef CGAL::Kd_tree<ST> Tree;
typedef CGAL::Fuzzy_sphere<ST> FSphere;
//typedef CGAL::Fuzzy_iso_box<ST> Fbox;

bool comp_sph(const Sphere_3 &a, const Sphere_3 &b) {
  return a.squared_radius()>b.squared_radius();
}

int main(int argc, char *argv[])
{
  std::string input=argv[1];
  std::string output=argv[2];
  FT lowRange=0;
  FT highRange=2592;
  FT cpyRange=100;
  FT outRange=50;
  FT boxsize=highRange-lowRange;
  FT radius;
  FT Vol;

  std::vector<Point> P;
  std::vector<Point> Pin;
  std::vector<Sphere_3> SP;
  std::map<Point_3,Sphere_3> index;
  std::map<Point_3,Sphere_3>::iterator it;
  Point_3 PP[4];
  Point_3 Ptmp;
  Sphere_3 Sph, Sph1;
  std::size_t size,i,k;
  long j;

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
  oFile.precision(10);

  for(Delaunay::Finite_cells_iterator cell=T.finite_cells_begin();cell!=T.finite_cells_end();cell++) {
    for(i=0;i<4;i++)
      PP[i] = cell->vertex(i)->point();
    Sph = Sphere_3(PP[0],PP[1],PP[2],PP[3]);
    Ptmp=Sph.center();
    if(Ptmp.x()<lowRange-outRange || Ptmp.x()>=highRange+outRange ||
        Ptmp.y()<lowRange-outRange || Ptmp.y()>=highRange+outRange ||
        Ptmp.z()<lowRange-outRange || Ptmp.z()>=highRange+outRange)
      continue;
    if(Ptmp.x()<lowRange || Ptmp.x()>=highRange ||
        Ptmp.y()<lowRange || Ptmp.y()>=highRange ||
        Ptmp.z()<lowRange || Ptmp.z()>=highRange)
      Vol=0;
    else
      Vol=CGAL::abs(CGAL::volume(PP[0],PP[1],PP[2],PP[3]));
    oFile<<Sph.center()<<" "<<CGAL::sqrt(CGAL::to_double(Sph.squared_radius()))<<" "<<Vol<<std::endl;
//    boxsize=CGAL::sqrt(CGAL::to_double(Sph.squared_radius()));
//    R.push_back(boxsize);
/*    SP.push_back(Sph);
    P.push_back(PP[0]);
    index[PP[0]]=Sph;
*/
  }
/*
  Tree tree(P.begin(),P.end());
  std::sort(SP.begin(),SP.end(),comp_sph);

  for(i=0;i<SP.size();i++) {
    Sph=SP[i];
    radius=CGAL::sqrt(CGAL::to_double(Sph.squared_radius()));
    oFile<<Sph.center()<<" "<<CGAL::make_sqrt(Sph.squared_radius())<<std::endl;

    FSphere fs(Sph.center(),radius*2,0);
    Pin.clear();
    tree.search(back_inserter(Pin),fs);
    //std::cout<<Pin.size()<<std::endl;;
    for(k=0;k<Pin.size();k++) {
      it=index.find(Pin[k]);
      if(it==index.end()) continue;
      Sph1=it->second;
      boxsize=CGAL::sqrt(CGAL::to_double(Sph1.squared_radius()));
      if(CGAL::squared_distance(Sph.center(),Pin[k])<CGAL::square(radius+boxsize))
        SP.erase(std::remove(SP.begin(),SP.end(),Sph1),SP.end());
    }
  }
*/

  oFile.close();

  std::cout<<"Finished successfully!"<<std::endl;
  return 0;
}

