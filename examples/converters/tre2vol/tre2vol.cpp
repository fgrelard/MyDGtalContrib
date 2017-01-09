/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file vol2sdp.cpp
 * @ingroup conerters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/07/21
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "shapes/Ball.h"
using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  std::vector<int> domainCoords;
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "Sequence of 3d Discrete points (.sdp) " )
    ("output,o", po::value<std::string>(), "Vol file  (.vol, .longvol, .pgm3d) " )
    ("foregroundVal,f", po::value<int>()->default_value(128), "value which will represent the foreground object in the resulting image (default 128)")
    ("invertY", "Invert the Y axis (image flip in the y direction)")
    ("backgroundVal,b", po::value<int>()->default_value(0), "value which will represent the background outside the  object in the resulting image (default 0)")
	  ("spacingX,x", po::value<double>()->default_value(1), "value which will represent the foreground object in the resulting image (default 128)")
	  ("spacingY,y", po::value<double>()->default_value(1), "value which will represent the foreground object in the resulting image (default 128)")
	  ("spacingZ,z", po::value<double>()->default_value(1), "value which will represent the foreground object in the resulting image (default 128)")
    ("domain,d",  po::value<std::vector <int> >(&domainCoords)->multitoken(), "The domain of the resulting image xmin ymin zmin xmax ymax zmax ");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
	  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [input-file] [output]\n"
		<< "Convert volumetric  file into a digital set of points from a given threshold."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.p3d \n";
      return 0;
    }
  if(! vm.count("input") ||! vm.count("output") )
    {
      trace.error() << " Input/ output filename and domain are needed to be defined" << endl;      
      return 0;
    }
  
  Z3i::Point ptLower(domainCoords[0],domainCoords[1], domainCoords[2]);
  Z3i::Point ptUpper(domainCoords[3],domainCoords[4], domainCoords[5]);
  Image3D::Domain imageDomain(ptLower, ptUpper);
  
  string inputSDP = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  int foregroundVal = vm["foregroundVal"].as<int>();
  int backgroundVal = vm["backgroundVal"].as<int>();
  double sx = vm["spacingX"].as<double>();
  double sy = vm["spacingY"].as<double>();
  double sz = vm["spacingZ"].as<double>();
  
  vector<unsigned int> vPos;
  vPos.push_back(0);
  vPos.push_back(1);
  vPos.push_back(2);
  trace.info() << "Reading input SDP file: " << inputSDP  ;
  typedef PointVector<4, double> TrePoint;
  std::vector<TrePoint> vectPoints=  PointListReader<TrePoint>::getPointsFromFile(inputSDP); 
  trace.info() << " [done] " << std::endl ; 

  Image3D imageResult(imageDomain); 
  for(Image3D::Domain::ConstIterator iter = imageResult.domain().begin(); iter!= imageResult.domain().end();
      iter++){
    imageResult.setValue(*iter, backgroundVal);
  }    
  
  for(unsigned int i=0; i<vectPoints.size(); i++){
	  Z3i::Point p(vectPoints[i][0]*sx, vectPoints[i][1]*sy, vectPoints[i][2]*sz);
	  Ball<Z3i::Point> ball(p, vectPoints[i][3]);
	  if(vm.count("invertY")){
		  p[1]=ptUpper[1]-p[1];
	  }
	  if(imageResult.domain().isInside(p)){
		  imageResult.setValue(p, foregroundVal);
		  vector<Z3i::Point> pointsInBall = ball.pointsInBall();
		  for (auto it = pointsInBall.begin(), ite = pointsInBall.end(); it != ite; ++it) {
			  if(imageResult.domain().isInside(*it))
				  imageResult.setValue(*it, foregroundVal);
		  }
	  }else{
		  trace.warning() << "point " << p << " outside the domain (ignored in the resulting volumic image)" << std::endl;  
	  }
  }
  trace.info()<< "Exporting resulting volumic image ... ";
  GenericWriter<Image3D>::exportFile(outputFilename, imageResult);
  trace.info() << " [done]"<<std::endl;
  return 0;
  
}




