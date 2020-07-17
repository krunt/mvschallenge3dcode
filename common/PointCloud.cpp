#include "PointCloud.h"
#include <ogr_spatialref.h>
#include "CorrImage.h"

//#include <liblas/liblas.hpp>
//#include <liblas/point.hpp>
//#include <liblas/reader.hpp>
//#include <liblas/writer.hpp>
//#include <liblas/color.hpp>


static bool IsNorthernHemisphere(double Lat) 
{
   return Lat >= 0;
}


// Get the UTM zone for a given Lat/Long
static int GetUTMZone(double Lat, double Long) 
{
   double LongTemp;
   int    ZoneNumber;

   /* Make sure the longitude is between -180.00 .. 179.9 */
   LongTemp = fmod( Long, 360 );
   if ( LongTemp < -180 ) LongTemp += 360;
   else if ( LongTemp > 180 ) LongTemp -= 360;

   ZoneNumber = (int)(( LongTemp + 180 ) / 6) + 1;

   if ( Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 ) ZoneNumber = 32;

   /* Special zones for Svalbard */
   if ( Lat >= 72.0 && Lat < 84.0 )
   {
      if (      LongTemp >=  0.0 && LongTemp <  9.0 ) ZoneNumber = 31;
      else if ( LongTemp >=  9.0 && LongTemp < 21.0 ) ZoneNumber = 33;
      else if ( LongTemp >= 21.0 && LongTemp < 33.0 ) ZoneNumber = 35;
      else if ( LongTemp >= 33.0 && LongTemp < 42.0 ) ZoneNumber = 37;
   }

   return ZoneNumber;
}


void MyPointCloud::ExportToXYZ(const std::string & outName, const MyCorrImage & refImage)
{
   OGRSpatialReference pj_utm, pj_latlong;
   pj_latlong.importFromProj4("+proj=latlong +datum=WGS84");
   pj_utm.importFromProj4("+proj=utm +datum=WGS84");
   int zone = GetUTMZone(m_bounds.GetCenter().y, m_bounds.GetCenter().x);
   pj_utm.SetUTM(zone, IsNorthernHemisphere(m_bounds.GetCenter().y));

   OGRCoordinateTransformation* coordinate_transform;
   coordinate_transform = OGRCreateCoordinateTransformation(&pj_latlong, &pj_utm);

   std::fstream fs(outName, std::ios_base::out);
   for (int i = 0; i < m_image.rows; ++i) {
      for (int j = 0; j < m_image.cols; ++j) {
         double lon = m_bounds.m_minLon + j * m_resolutionX;
         double lat = m_bounds.m_minLat + i * m_resolutionY;
         double h = m_image.at<double>(i, j);

         // skip nans
         if (h != h) {
            continue;
         }

         double intensity = std::max(refImage.GetIntensityByLonLat(MyPoint2d(lon,lat),h),0.0);

         coordinate_transform->Transform(1, &lon, &lat);
         //fs << lon << " " << lat << " " << h << " " << intensity << "\n";
         char buf[1024];
         int sz = sprintf(buf, "%.3f %.3f %.3f %.3f\n", lon, lat, h, intensity);
         fs.write(buf,sz);
      }
   }
}


void MyPointCloud::ExportToLAS(const std::string & outName, const MyCorrImage & refImage)
{
   //std::ofstream ofs;
   //if (!liblas::Create(ofs, outName)) {
   //   throw std::runtime_error(std::string("Can not create ") + outName);
   //}

   //liblas::HeaderPtr header(new liblas::Header());

   //// Checks if the filename extension ends in z and enables compression if
   //// it does
   //if (outName.back() == 'z') {
   //   header->SetCompressed(true);
   //} else {
   //   header->SetCompressed(false); 
   //}

   //OGRSpatialReference pj_utm, pj_latlong;
   //pj_latlong.importFromProj4("+proj=latlong +datum=WGS84");
   //pj_utm.importFromProj4("+proj=utm +datum=WGS84");
   //int zone = GetUTMZone(m_bounds.GetCenter().y, m_bounds.GetCenter().x);
   //pj_utm.SetUTM(zone, IsNorthernHemisphere(m_bounds.GetCenter().y));

   //OGRCoordinateTransformation* coordinate_transform;
   //coordinate_transform = OGRCreateCoordinateTransformation(&pj_latlong, &pj_utm);

   //int cnt = m_image.rows*m_image.cols;
   //double min_x=1e30, min_y=1e30, min_z=1e30;
   //double max_x=-1e30, max_y=-1e30, max_z=-1e30;
   //for (int i = 0; i < m_image.rows; ++i) {
   //   for (int j = 0; j < m_image.cols; ++j) {
   //      double lon = m_bounds.m_minLon + j * m_resolutionX;
   //      double lat = m_bounds.m_minLat + i * m_resolutionY;
   //      double h = m_image.at<double>(i, j);
   //      coordinate_transform->Transform(1, &lon, &lat);
   //      double intensity = std::max(refImage.GetIntensityByLonLat(MyPoint2d(lon,lat),h),0.0);

   //      min_x = std::min<double>(min_x, lon);
   //      min_y = std::min<double>(min_y, lat);
   //      min_z = std::min<double>(min_z, h);

   //      max_x = std::max<double>(max_x, lon);
   //      max_y = std::max<double>(max_y, lat);
   //      max_z = std::max<double>(max_z, h);
   //   }
   //}


   //header->SetMax(max_x, max_y, max_z);
   //header->SetMin(min_x, min_y, min_z);


   //// Calculates an good scale and offset to reduce quatization error when
   //// the points are packed into the file (LAS format stores positions as
   //// 32 bit integers)
   //double scale_x = ((max_x - min_x)/2.0) / std::numeric_limits<int32_t>::max();
   //double scale_y = ((max_y - min_y)/2.0) / std::numeric_limits<int32_t>::max();
   //double scale_z = ((max_z - min_z)/2.0) / std::numeric_limits<int32_t>::max();

   //// Offset at the center minimizes quatization error
   //double offset_x = (max_x - min_x)/2.0 + min_x;
   //double offset_y = (max_y - min_y)/2.0 + min_y;
   //double offset_z = (max_z - min_z)/2.0 + min_z;

   //header->SetScale(scale_x, scale_y, scale_z);
   //header->SetOffset(offset_x, offset_y, offset_z);

   //// Specify the LAS format as 1.2
   //header->SetVersionMajor(1);
   //header->SetVersionMinor(2);

   //header->SetDataFormatId(liblas::ePointFormat0);
   //header->SetPointRecordsCount(cnt);

   //// Both these lines are needed to update the georefence from the spatial_ref
   //// stored in the point cloud
   //char* proj4_str;
   //pj_utm.exportToProj4(&proj4_str);
   //liblas::SpatialReference lasReference;
   //lasReference.SetProj4(proj4_str);
   //
   //header->SetSRS(const_cast<liblas::SpatialReference&>(lasReference));
   //header->SetGeoreference();

   //liblas::Writer writer(ofs, *header);

   //for (int i = 0; i < m_image.rows; ++i) {
   //   for (int j = 0; j < m_image.cols; ++j) {
   //      double lon = m_bounds.m_minLon + j * m_resolutionX;
   //      double lat = m_bounds.m_minLat + i * m_resolutionY;
   //      double h = m_image.at<double>(i, j);
   //      coordinate_transform->Transform(1, &lon, &lat);
   //      double intensity = std::max(refImage.GetIntensityByLonLat(MyPoint2d(lon,lat),h),0.0);

   //      liblas::Point point(header);
   //      point.SetCoordinates(lon,lat,h);
   //      point.SetIntensity(intensity);

   //      point.SetReturnNumber(0);
   //      point.SetNumberOfReturns(0);
   //      point.SetFlightLineEdge(0);

   //      writer.WritePoint(point);
   //   }
   //}

   //writer.WriteHeader();
}