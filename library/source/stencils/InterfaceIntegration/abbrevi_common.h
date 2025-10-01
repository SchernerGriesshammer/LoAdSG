/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

//////////////////////////////////////////////////////////
//                  abbbreviations
//////////////////////////////////////////////////////////

#ifndef ABBREVICOMMON_H_
#define ABBREVICOMMON_H_

#ifdef _MPI_PARALLEL
#include <mpi.h>
#endif 

#include <iostream>
#include <fstream>
using std::ofstream;
using std::cout;
using std::endl;

#include <complex>
#include <cmath>

#include <map>
#include <cassert>



#ifndef M_PI
#define M_PI 3.1415927
#endif

enum elementTyp { pointEl, edgeEl,  rectangleEl, triangleEl, quadrangleEl, hexahedronEl };


// 1D
//-------------
enum dir1D { Ldir1D,  Rdir1D };
enum Ort1D { LOrt, MOrt, ROrt };

// 2D
//-------------
// directions in 2D
enum dir2D { Wdir2D, Edir2D, Sdir2D, Ndir2D };
// direction of sons in 2D
enum dir2D_sons { SWdir2D, SEdir2D, NWdir2D, NEdir2D };
// local directions of hexahedra with fixed SD-direction
enum local_dir_of_hex2D{ ND, ST, NT };

// 3D
//-------------
// directions in 3D
enum dir3D { Wdir3D, Edir3D, Sdir3D, Ndir3D, Ddir3D, Tdir3D };
// directions fuer Soehne, oder Ecken einer Zellen
enum dir3D_sons { WSDdir3D, ESDdir3D, WNDdir3D, ENDdir3D, WSTdir3D, ESTdir3D, WNTdir3D, ENTdir3D };
// local directions of hexahedra with fixed SWD-direction
enum local_dir_of_hex3D{ SED, NWD, NED, SWT, SET, NWT, NET };
// main_directions in 3D
enum main_dir_3D { WE_dir, SN_dir, DT_dir };

// edges of a cell
enum Edges_cell { SDed, NDed, STed, NTed, WDed, EDed, WTed, ETed,
                  SWed, SEed, NWed, NEed };

// gives back eta or xi or phi		  
inline double find_p ( Edges_cell ed, double eta, double xi,double phi ) {
  if ( ed<=NTed ) return eta;
  if ( ed<=ETed ) return xi;
  if ( ed>NEed ) cout << " error in find_p!" << endl;
  return phi;
}


// corners of an edge of a cell 
//------------------------------
inline dir3D_sons Transform(Edges_cell ec, dir1D d) {
  enum dir3D dir;
  dir = (dir3D) ((ec >> 2) << 1);
  if(dir<Sdir3D) return (dir3D_sons)(d+4*((ec>>1)&1)+2*(ec&1));
  if(dir>Ndir3D) return (dir3D_sons)(2*((ec>>1)&1)+(ec&1)+4*d);
  return (dir3D_sons)(4*((ec>>1)&1)+2*d+(ec&1));
}

#endif // ABBREVICOMMON_H_
