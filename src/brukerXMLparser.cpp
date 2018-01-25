/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2018 Pere Rafols Soler
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#include "pugixml.hpp"
#include <Rcpp.h>
#include <string>
#include <stdlib.h>   //for the atoi function

using namespace Rcpp;
using namespace pugi;

//' ParseBrukerXML.
//'
//' Reads a Bruker's xml file exported using fleximaging.
//' A list is returned where each element in the list is named according the ROI name.
//' Each element in the list consists in a data.frame with the pixels XY coordinates inside each ROI.
//'
//' @param xml_path the full path where XML file is stored.
//'
//' @return ROI pixel coordinates arranged in a named list.
//' 
// [[Rcpp::export]]
List CparseBrukerXML( String xml_path )
{
  std::string sName;
  std::string sValue;
  xml_document doc;
  xml_parse_result  result = doc.load_file(xml_path.get_cstring());
  if (!result)
  {
    Rcout << "XML [" << xml_path.get_cstring() << "] parsed with errors, attr value: [" << doc.child("node").attribute("attr").value() << "]" << std::endl;
    Rcout << "Error description: " << result.description() << std::endl;
    Rcout << "Error offset: " << result.offset << " (error at [..." << (xml_path.get_cstring() + result.offset) << "]" << std::endl;
    return( List::create(Named("Error") = "XML parse error") );
  }
  
  //Load the root node
  xml_node root = doc.child("ClinProtSpectraImport");
  if( root == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no ClinProtSpectraImport node found") );
  }
  
  const int iNumRois = root.select_nodes("Class").size();
  StringVector sRoiNames = StringVector(iNumRois);
  List lstRois = List( iNumRois );
  lstRois.names() = sRoiNames; //link list names reference to sRoiNames
  
  //Load every ROI node (Named Class in the XML)
  int iRoi = 0;
  for (xml_node ROI = root.child("Class"); ROI; ROI = ROI.next_sibling("Class"))
  {
    if( iRoi >= iNumRois )
    {
      return( List::create(Named("Error") = "XML parse error: ROI counter overflow") );
    }
    
    if( ROI == NULL )
    {
      return( List::create(Named("Error") = "XML parse error: no Class node found") );
    }
    
    //Get ROI name
    sName = ROI.attribute("Name").value();
    sRoiNames[iRoi] = sName;
    
    //Parse each pixel in current ROI
    const int iNumPixelsInRoi = ROI.select_nodes("Element").size();
    int iPixel = 0;
    std::size_t xPos, yPos;
    IntegerVector iX(iNumPixelsInRoi);
    IntegerVector iY(iNumPixelsInRoi);
    StringVector sSpots(iNumPixelsInRoi);
    for (xml_node pixel = ROI.child("Element"); pixel; pixel = pixel.next_sibling("Element"))
    {
      if( iPixel >= iNumPixelsInRoi )
      {
        return( List::create(Named("Error") = "XML parse error: pixel counter overflow") );
      }
      
      if( pixel == NULL )
      {
        return( List::create(Named("Error") = "XML parse error: no Element node found") );
      }
      
      sValue = pixel.attribute("Spot").value();
      sSpots[iPixel] = sValue;
      xPos = sValue.find('X');
      yPos = sValue.find('Y');
      if( xPos == std::string::npos || yPos == std::string::npos ||
          xPos >= yPos || yPos >= (sValue.length() - 1))
      {
        return( List::create(Named("Error") = "XML parse error: no valid Spot attribute node found") );
      }
      
      //Parse pixel coordinates
      iX[iPixel] = atoi( sValue.substr(xPos + 1, yPos - xPos - 1).c_str() );
      iY[iPixel] = atoi( sValue.substr(yPos + 1, sValue.length() - yPos - 1).c_str() );
      iPixel++;
    }
    //Build the roi list
    lstRois[iRoi] =  List::create( Named("name") = sName.c_str(),
                                   Named("pos") = DataFrame::create( Named("x") = iX, Named("y") = iY),
                                   Named("spots") = sSpots );
      
    
    iRoi++;
  }
  
  return lstRois;
}

