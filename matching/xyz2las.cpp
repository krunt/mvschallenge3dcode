#include <liblas/liblas.hpp>
#include <liblas/point.hpp>
#include <liblas/reader.hpp>
#include <liblas/writer.hpp>
#include <liblas/color.hpp>

int main(int argc, char **argv)
{
   if (argc != 3) {
      return 1;
   }

   std::string inName = argv[1];
   std::string outName = argv[2];

   std::ofstream ofs;
   if (!liblas::Create(ofs, outName)) {
      throw std::runtime_error(std::string("Can not create ") + outName);
   }

   liblas::HeaderPtr header(new liblas::Header());

   // Checks if the filename extension ends in z and enables compression if
   // it does
   if (outName.back() == 'z') {
      header->SetCompressed(true);
   } else {
      header->SetCompressed(false); 
   }

   OGRSpatialReference pj_utm;
   pj_utm.importFromProj4("+proj=utm +datum=WGS84");
   pj_utm.SetUTM(21, false);

   double min_x=1e30, min_y=1e30, min_z=1e30;
   double max_x=-1e30, max_y=-1e30, max_z=-1e30;

   {
      int cnt = 0;
      std::fstream ifs(inName, std::ios_base::in);
      for (;;) {
         double x, y, z, intensity;
         ifs >> x >> y >> z >> intensity;
         if (ifs.eof()) {
            break;
         }

         min_x = std::min<double>(min_x, lon);
         min_y = std::min<double>(min_y, lat);
         min_z = std::min<double>(min_z, h);

         max_x = std::max<double>(max_x, lon);
         max_y = std::max<double>(max_y, lat);
         max_z = std::max<double>(max_z, h);

         ++cnt;
      }
   }

   header->SetMax(max_x, max_y, max_z);
   header->SetMin(min_x, min_y, min_z);


   // Calculates an good scale and offset to reduce quatization error when
   // the points are packed into the file (LAS format stores positions as
   // 32 bit integers)
   double scale_x = ((max_x - min_x)/2.0) / std::numeric_limits<int32_t>::max();
   double scale_y = ((max_y - min_y)/2.0) / std::numeric_limits<int32_t>::max();
   double scale_z = ((max_z - min_z)/2.0) / std::numeric_limits<int32_t>::max();

   // Offset at the center minimizes quatization error
   double offset_x = (max_x - min_x)/2.0 + min_x;
   double offset_y = (max_y - min_y)/2.0 + min_y;
   double offset_z = (max_z - min_z)/2.0 + min_z;

   header->SetScale(scale_x, scale_y, scale_z);
   header->SetOffset(offset_x, offset_y, offset_z);

   // Specify the LAS format as 1.2
   header->SetVersionMajor(1);
   header->SetVersionMinor(2);

   header->SetDataFormatId(liblas::ePointFormat0);
   header->SetPointRecordsCount(cnt);

   // Both these lines are needed to update the georefence from the spatial_ref
   // stored in the point cloud
   char* proj4_str;
   pj_utm.exportToProj4(&proj4_str);
   liblas::SpatialReference lasReference;
   lasReference.SetProj4(proj4_str);
   
   header->SetSRS(const_cast<liblas::SpatialReference&>(lasReference));
   header->SetGeoreference();

   liblas::Writer writer(ofs, *header);

   std::fstream ifs(inName, std::ios_base::in);
   for (int i = 0; i < cnt; ++i) {
      double x, y, z, intensity;
      ifs >> x >> y >> z >> intensity;
      if (ifs.eof()) {
         break;
      }

      liblas::Point point(header);
      point.SetCoordinates(x,y,z);
      point.SetIntensity(intensity);

      point.SetReturnNumber(0);
      point.SetNumberOfReturns(0);
      point.SetFlightLineEdge(0);

      writer.WritePoint(point);
   }

   writer.WriteHeader();
   return 0;
}