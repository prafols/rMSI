/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2019 Pere Rafols Soler
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

#include "common_methods.h"

std::string parse_xml_uuid(std::string uuid)
{
  std::size_t ipos = uuid.find('{');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  do
  {
    ipos = uuid.find('-'); 
    if( ipos != std::string::npos )
    {
      uuid.erase(ipos, 1);
    } 
  } while ( ipos != std::string::npos);
  ipos = uuid.find('}');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  for( unsigned int i=0; i < uuid.length(); i++)
  {
    uuid[i] = toupper(uuid[i]);
  }
  return uuid;
}
