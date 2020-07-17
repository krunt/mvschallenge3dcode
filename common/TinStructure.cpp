
#include  "TinStructure.h"

#ifdef DEBUG_MATCHER
#pragma optimize("",off)
#endif

void TinStructure::Construct(const std::vector<cv::Point3d> & pnts, const std::vector<std::pair<int,int>> & edges)
{
   m_pnts = pnts;
   m_edges = edges;

   struct triangulateio in;
   struct triangulateio &out = m_structure;
   
   memset(&in,0,sizeof(in));
   memset(&out,0,sizeof(out));
   
   in.numberofpoints = (int)pnts.size();
   in.pointlist = (REAL*)malloc(in.numberofpoints * 2 * sizeof(REAL));
   in.numberofpointattributes = 1;
   in.pointattributelist = (REAL*)malloc(in.numberofpoints * sizeof(REAL));
   
   for (int i = 0; i < in.numberofpoints; ++i) {
      in.pointlist[i*2+0] = pnts[i].x;
      in.pointlist[i*2+1] = pnts[i].y;
      in.pointattributelist[i] = pnts[i].z;
   }
   
   in.pointmarkerlist = NULL;
   in.holelist = NULL;
   in.numberofholes = 0;
   
   in.numberofsegments = (int)edges.size();
   in.segmentlist = (int*)malloc(edges.size()*2*sizeof(int));
   in.segmentmarkerlist  = (int*)malloc(edges.size()*sizeof(int));
   
   for (int i = 0; i < edges.size(); ++i) {
      in.segmentlist[2*i+0] = edges[i].first;
      in.segmentlist[2*i+1] = edges[i].second;
      in.segmentmarkerlist[i] = i+2;
   }
   
   out.pointlist = NULL;
   out.pointattributelist = NULL;
   out.pointmarkerlist = NULL;
   out.trianglelist = NULL;
   out.neighborlist = NULL;
   out.segmentlist = NULL;
   out.segmentmarkerlist = NULL;
   out.edgelist = NULL;
   out.edgemarkerlist = NULL;

   m_meshStruct = NULL;
   m_behavStruct = NULL;
   
   triangulate("necpz", &in, &out, NULL, &m_meshStruct, &m_behavStruct);
}

// lon/lat
bool TinStructure::GetTriangleByPoint(const cv::Point2d & pnt, std::vector<TinPoint> & triangle)
{
   REAL queryPoint[2];

   if (!m_meshStruct || !m_behavStruct) {
      return false;
   }

   queryPoint[0] = pnt.x;
   queryPoint[1] = pnt.y;

   triangle.clear();

   int outindices[6];
   if (delaunaylocate(m_meshStruct, m_behavStruct, queryPoint, outindices) == 0) {
      for (int i = 0; i < 3; ++i) {
         REAL * pt = &m_structure.pointlist[2*outindices[i]];
         REAL altitude = m_structure.pointattributelist[outindices[i]];
         bool isLieOnSegment = m_structure.pointmarkerlist[outindices[i]] > 1;
         triangle.push_back(TinPoint(cv::Point3d(pt[0],pt[1],altitude),isLieOnSegment));
      }

      return true;
   }

   return false;
}


static void FetchPointNeighbors(int pointIdx, int triIdx, struct triangulateio & inp, std::set<int> & visited, std::set<int> & neigh)
{
   if (triIdx==-1) {
      return;
   }

   if (visited.count(triIdx) > 0) {
      return;
   }

   int p0 = inp.trianglelist[triIdx*inp.numberofcorners+0];
   int p1 = inp.trianglelist[triIdx*inp.numberofcorners+1];
   int p2 = inp.trianglelist[triIdx*inp.numberofcorners+2];

   bool needRecurse = false;
   if (p0==pointIdx) { 
      needRecurse = true; neigh.insert(triIdx);
   }
   if (p1==pointIdx) { 
      needRecurse = true; neigh.insert(triIdx);
   }
   if (p2==pointIdx) { 
      needRecurse = true; neigh.insert(triIdx);
   }

   visited.insert(triIdx);

   if (needRecurse) {
      for (int i = 0; i < 3; ++i) {
         if (inp.neighborlist[triIdx*3+i] != -1) {
            FetchPointNeighbors(pointIdx, inp.neighborlist[triIdx*3+i], inp, visited, neigh);
         }
      }
   }
}

//#pragma optimize("",off)
bool TinStructure::GetNeighborsForPoint(const cv::Point2d & pnt, std::vector<TinPoint> & neighbors)
{
   REAL queryPoint[2];

   if (!m_meshStruct || !m_behavStruct) {
      return false;
   }

   queryPoint[0] = pnt.x;
   queryPoint[1] = pnt.y;

   int outindices[6];
   if (delaunaylocate(m_meshStruct, m_behavStruct, queryPoint, outindices) == 0) {

      int pointIdx = -1;
      for (int i = 0; i < 3; ++i) {
         REAL * pt = &m_structure.pointlist[2*outindices[i]];
         if (abs(pt[0]-pnt.x)< MYEPS && abs(pt[1]-pnt.y)< MYEPS) {
            pointIdx = outindices[i];
         }
      }

      if (pointIdx == -1) {
         return false;
      }

      std::set<int> neigh, visited, points;
      FetchPointNeighbors(pointIdx, outindices[3], m_structure, visited, neigh);
      FetchPointNeighbors(pointIdx, outindices[4], m_structure, visited, neigh);
      FetchPointNeighbors(pointIdx, outindices[5], m_structure, visited, neigh);

      for (auto triIndex: neigh) {
         int p0 = m_structure.trianglelist[triIndex*m_structure.numberofcorners+0];
         int p1 = m_structure.trianglelist[triIndex*m_structure.numberofcorners+1];
         int p2 = m_structure.trianglelist[triIndex*m_structure.numberofcorners+2];
         int pts[3] = {p0,p1,p2};
         for (int j = 0; j < 3; ++j) {
            if (pts[j] != pointIdx) {
               points.insert(pts[j]);
            }
         }
      }

      for (auto pointIdx: points) {
         REAL * pt = &m_structure.pointlist[2*pointIdx];
         REAL altitude = m_structure.pointattributelist[pointIdx];
         bool isLieOnSegment = m_structure.pointmarkerlist[pointIdx] > 1;
         neighbors.push_back(TinPoint(cv::Point3d(pt[0],pt[1],altitude),isLieOnSegment));
      }

      return true;
   }

   return false;
}
//#pragma optimize("",on)


static void CalcBarycentric(cv::Point2d p, const std::vector<TinPoint> & abc, float &u, float &v, float &w)
{
   const cv::Point3d &a = abc[0].m_pnt;
   const cv::Point3d &b = abc[1].m_pnt;
   const cv::Point3d &c = abc[2].m_pnt;

   cv::Point2f v0, v1, v2;
   v0 = cv::Point2f(float(b.x-a.x),float(b.y-a.y));
   v1 = cv::Point2f(float(c.x-a.x),float(c.y-a.y));
   v2 = cv::Point2f(float(p.x-a.x),float(p.y-a.y));
   
   float d00 = v0.dot(v0);
   float d01 = v0.dot(v1);
   float d11 = v1.dot(v1);
   float d20 = v2.dot(v0);
   float d21 = v2.dot(v1);
   float denom = d00 * d11 - d01 * d01;
   v = (d11 * d20 - d01 * d21) / denom;
   w = (d00 * d21 - d01 * d20) / denom;
   u = 1.0f - v - w;
}


#ifdef DEBUG_MATCHER
#pragma optimize("",off)
#endif
void TinStructure::Interpolate(MyPointCloud & out)
{
   int ww = out.GetWidth();
   int hh = out.GetHeight();

   std::vector<TinPoint> triangle;
   for (int i = 0; i < hh; ++i) {
      for (int j = 0; j < ww; ++j) {
         cv::Point2d pnt = cv::Point2d(j,i);
         if (GetTriangleByPoint(pnt, triangle)) {

            float u, v, w;
            CalcBarycentric(pnt, triangle, u, v, w);

            float h = u * (float)triangle[0].m_pnt.z 
               + v * (float)triangle[1].m_pnt.z + w * (float)triangle[2].m_pnt.z;
            out.SetAltitude(i, j, (double)h);
         }
      }
   }
}


void TinStructure::DebugExport(MyPointCloud & out)
{
   int ww = out.GetWidth();
   int hh = out.GetHeight();

   double kScale=1.0;
   cv::Mat img = cv::Mat::zeros(hh/kScale, ww/kScale, CV_8UC1);

   for (int i = 0; i < m_structure.numberoftriangles; ++i) {
      for (int j = 1; j <= 3; ++j) {
         int from = j - 1;
         int to = j % 3;

         from = m_structure.trianglelist[i*3+from];
         to = m_structure.trianglelist[i*3+to];

         cv::Point fromPoint = cv::Point(m_structure.pointlist[from*2+0]/kScale,m_structure.pointlist[from*2+1]/kScale);
         cv::Point toPoint = cv::Point(m_structure.pointlist[to*2+0]/kScale,m_structure.pointlist[to*2+1]/kScale);

         cv::line(img, fromPoint, toPoint, cv::Scalar(64,255,255));
      }
   }

   for (int i = 0; i < m_structure.numberofsegments; ++i) {
      int from = m_structure.segmentlist[2*i+0];
      int to = m_structure.segmentlist[2*i+1];

      cv::Point fromPoint = cv::Point(m_structure.pointlist[from*2+0]/kScale,m_structure.pointlist[from*2+1]/kScale);
      cv::Point toPoint = cv::Point(m_structure.pointlist[to*2+0]/kScale,m_structure.pointlist[to*2+1]/kScale);

      cv::line(img, fromPoint, toPoint, cv::Scalar(128,255,255));
   }

   for (int i = 0; i < m_structure.numberofpoints; ++i) {
      cv::circle(img, cv::Point(m_structure.pointlist[i*2+0]/kScale,m_structure.pointlist[i*2+1]/kScale), 1, cv::Scalar(255,255,255));
   }

   for (int i = 0; i < hh; ++i) {
      for (int j = 0; j < ww; ++j) {
         uchar v = img.at<uchar>(i,j);
         if (v>0) {
            out.SetAltitude(i, j, 40.0);
         }
      }
   }
}


void TinStructure::Visualize(int width, int height) const
{
   cv::namedWindow("myWin");
   
   double kScale=2.0;
   cv::Mat img = cv::Mat::zeros(height/kScale, width/kScale, CV_8UC3);
   for (int i = 0; i < m_structure.numberofpoints; ++i) {
      cv::circle(img, cv::Point(m_structure.pointlist[i*2+0]/kScale,m_structure.pointlist[i*2+1]/kScale), 2, cv::Scalar(255,255,255));
   }
   for (int i = 0; i < m_structure.numberoftriangles; ++i) {
      for (int j = 1; j <= 3; ++j) {
         int from = j - 1;
         int to = j % 3;
   
         from = m_structure.trianglelist[i*3+from];
         to = m_structure.trianglelist[i*3+to];
   
         cv::Point fromPoint = cv::Point(m_structure.pointlist[from*2+0]/kScale,m_structure.pointlist[from*2+1]/kScale);
         cv::Point toPoint = cv::Point(m_structure.pointlist[to*2+0]/kScale,m_structure.pointlist[to*2+1]/kScale);
   
         cv::line(img, fromPoint, toPoint, cv::Scalar(255,255,255));
      }
   }
   
   for (int i = 0; i < m_structure.numberofsegments; ++i) {
      int from = m_structure.segmentlist[2*i+0];
      int to = m_structure.segmentlist[2*i+1];
   
      cv::Point fromPoint = cv::Point(m_structure.pointlist[from*2+0]/kScale,m_structure.pointlist[from*2+1]/kScale);
      cv::Point toPoint = cv::Point(m_structure.pointlist[to*2+0]/kScale,m_structure.pointlist[to*2+1]/kScale);
   
      cv::line(img, fromPoint, toPoint, cv::Scalar(0,255,255));
   }
   cv::imshow("myWin", img);
   cv::waitKey();
}

#ifdef DEBUG_MATCHER
#pragma optimize("",on)
#endif