#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
// includes from get_triangle.cxx
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

#include <vtkDataSetWriter.h>
#include "shading.cxx"
#include <cmath>

#define TEST 0

using std::cerr;
using std::endl;

double ceil441(double f)
{
  return ceil(f-0.00001);
}

double floor441(double f)
{
  return floor(f+0.00001);
}

double lerp(double v, double v1, double v2, double A, double B) {

  if (B == A) {
    return A; 
  }
  if (v2 == v1) { 
    return B;
  }

  return A + ((v-v1)*(B-A)) / (v2 - v1);
}

vtkImageData *NewImage(int width, int height)
{
  vtkImageData *img = vtkImageData::New();
  img->SetDimensions(width, height, 1);
  img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

  return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
  std::string full_filename = filename;
  full_filename += ".png";
  vtkPNGWriter *writer = vtkPNGWriter::New();
  writer->SetInputData(img);
  writer->SetFileName(full_filename.c_str());
  writer->Write();
  writer->Delete();
}

class Color {
  public:
    double R;
    double G;
    double B;

    // Are these frowned upon?
    Color operator+(const Color& b) {
      Color color;
      color.R = this->R + b.R;
      color.G = this->G + b.G;
      color.B = this->B + b.B;

      return color;
    }

    Color operator*(double d) {
      Color color;
      color.R = this->R * d;
      color.G = this->G * d;
      color.B = this->B * d;

      return color;
    }

};

Color lerp_color(double v, double v1, double v2, Color a, Color b) {

  Color c;

  c.R = lerp(v, v1, v2, a.R, b.R);
  c.G = lerp(v, v1, v2, a.G, b.G);
  c.B = lerp(v, v1, v2, a.B, b.B);

  return c;
}

class Normal {
  public:
    double X;
    double Y;
    double Z;

    Normal normalize(void) 
    {
      Normal n;

      double length = sqrt( (X*X + Y*Y + Z*Z));
      n.X = X / length;
      n.Y = Y / length;
      n.Z = Z / length;

      return n;
    }
    Normal operator*(double d) {
      Normal n;
      n.X = this->X * d;
      n.Y = this->Y * d;
      n.Z = this->Z * d;

      return n;
    }
};

double dot_product(Normal a, Normal b) {

  return (a.X * b.X) + (a.Y * b.Y) + (a.Z * b.Z);
}


Normal lerp_normal(double v, double v1, double v2, Normal a, Normal b) {

  Normal n;

  n.X = lerp(v, v1, v2, a.X, b.X);
  n.Y = lerp(v, v1, v2, a.Y, b.Y);
  n.Z = lerp(v, v1, v2, a.Z, b.Z);

  return n.normalize();
}


class Line {
  public:
    Line(double x1, double y1, double x2, double y2);
    double getX(double y);
    void printLine();

  private:
    double m;
    double b;

    int vert;
    double vert_x;
};

Line::Line(double x1, double y1, double x2, double y2) {
  if (x1 == x2) {
    vert = 1;
    vert_x = x1;
    return;
  }
  else {
    vert = 0;
  }

  m = (y2 - y1) / (x2 - x1);
  b = y1 - m*x1;

  return;
}

double Line::getX(double y) { 
  if (vert) {
    return vert_x;
  }
  //else if (m == 0) {
  //  return 0;
  //}
  else {
    return (y - b) / m;
  }
}

void Line::printLine() { 
  if (vert) {
    printf("vertical line at x = %f\n", vert_x);
  }
  else {
    printf("y = %f * x + %f\n", m, b);
  }
}

class Triangle
{
  public:
    double  X[3];
    double  Y[3];
    double  Z[3];
    Color   colors[3];
    Normal  normals[3];

    void print_triangle();
};

void Triangle::print_triangle() {
  for (int i = 0; i < 3; i++) {
    printf("Tri %d: X: %f, Y: %f, Z: %f\n R: %f, G: %f, B: %f\n, Nx: %f, Ny: %f, Nz: %f\n", 
        i, X[i], Y[i], Z[i], colors[i].R, colors[i].G, colors[i].B, normals[i].X, normals[i].Y, normals[i].Z);
  }
}

double CalculatePhongShading(LightingParameters &, Normal light, Normal normal)
{
  double dp = dot_product(light.normalize(), normal.normalize());

  // Ambient switch
  double ambient = 1;

  // Calulate two-sided for diffuse
  double diffuse = std::abs(dp); 

  // Calculate one-sided for specular
  Normal R;  // technically not a normal...
  R.X = 2 * dp * normal.X - light.X;
  R.Y = 2 * dp * normal.Y - light.Y;
  R.Z = 2 * dp * normal.Z - light.Z;

  Normal V;
  V.X = 0;
  V.Y = 0;
  V.Z = -1;

  double specular = pow(dot_product(V, R), lp.alpha);
  specular = fmax(0, specular); 

  return lp.Ka*ambient + lp.Kd*diffuse + lp.Ks*specular;
}

class Screen
{
  public:
    unsigned char *buffer;
    double *Zgrid;
    int width, height;

    // would some methods for accessing and setting pixels be helpful?
    void ImageColor(int r, int c, double z, Color color, Normal normal);
};

void Screen::ImageColor(int r, int c, double z, Color color, Normal normal) {
  if (r >= height || r < 0 || c >= width || c < 0) return;

  if (z < Zgrid[r*width + c]) return;

  Zgrid[r*width + c] = z;

  Normal light;
  light.X = lp.lightDir[0];
  light.Y = lp.lightDir[1];
  light.Z = lp.lightDir[2];

  double shade = CalculatePhongShading(lp, light, normal);

  color = color * shade;
  color.R = fmin(color.R, 1);
  color.G = fmin(color.G, 1);
  color.B = fmin(color.B, 1);

  color.R = fmax(color.R, 0);
  color.G = fmax(color.G, 0);
  color.B = fmax(color.B, 0);

  int offset = 3*(r*width + c);
  buffer[offset]     = ceil441(color.R * 255);
  buffer[offset + 1] = ceil441(color.G * 255);
  buffer[offset + 2] = ceil441(color.B * 255);
}

#if TEST
std::vector<Triangle> GetTriangles(void)
{
  std::vector<Triangle> rv(1);

  rv[0].X[0] = 100;
  rv[0].X[1] = 100;
  rv[0].X[2] = 600;

  rv[0].Y[0] = 100;
  rv[0].Y[1] = 800;
  rv[0].Y[2] = 800;

  rv[0].colors[0].R = 1;
  rv[0].colors[0].G = 0;
  rv[0].colors[0].B = 0;

  rv[0].colors[1].R = 0;
  rv[0].colors[1].G = 1;
  rv[0].colors[1].B = 0;

  rv[0].colors[2].R = 0;
  rv[0].colors[2].G = 0;
  rv[0].colors[2].B = 1;

  return rv;
}

#else
  std::vector<Triangle>
GetTriangles(void)
{
  vtkPolyDataReader *rdr = vtkPolyDataReader::New();
  rdr->SetFileName("proj1e_geometry.vtk");
  cerr << "Reading" << endl;
  rdr->Update();
  cerr << "Done reading" << endl;
  if (rdr->GetOutput()->GetNumberOfCells() == 0)
  {
    cerr << "Unable to open file!!" << endl;
    exit(EXIT_FAILURE);
  }
  vtkPolyData *pd = rdr->GetOutput();

  int numTris = pd->GetNumberOfCells();
  vtkPoints *pts = pd->GetPoints();
  vtkCellArray *cells = pd->GetPolys();
  vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
  double *color_ptr = var->GetPointer(0);
  //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
  //float *color_ptr = var->GetPointer(0);
  vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
  float *normals = n->GetPointer(0);
  std::vector<Triangle> tris(numTris);
  vtkIdType npts;
  vtkIdType *ptIds;
  int idx;
  for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
  {
    if (npts != 3)
    {
      cerr << "Non-triangles!! ???" << endl;
      exit(EXIT_FAILURE);
    }
    double *pt = NULL;
    pt = pts->GetPoint(ptIds[0]);
    tris[idx].X[0] = (pt[0]+10)*50.0;
    tris[idx].Y[0] = (pt[1]+10)*50.0;
    tris[idx].Z[0] = (pt[2]-10)*0.05;
    tris[idx].normals[0].X = normals[3*ptIds[0]+0];
    tris[idx].normals[0].Y = normals[3*ptIds[0]+1];
    tris[idx].normals[0].Z = normals[3*ptIds[0]+2];
    pt = pts->GetPoint(ptIds[1]);
    tris[idx].X[1] = (pt[0]+10)*50.0;
    tris[idx].Y[1] = (pt[1]+10)*50.0;
    tris[idx].Z[1] = (pt[2]-10)*0.05;
    tris[idx].normals[1].X = normals[3*ptIds[1]+0];
    tris[idx].normals[1].Y = normals[3*ptIds[1]+1];
    tris[idx].normals[1].Z = normals[3*ptIds[1]+2];
    pt = pts->GetPoint(ptIds[2]);
    tris[idx].X[2] = (pt[0]+10)*50.0;
    tris[idx].Y[2] = (pt[1]+10)*50.0;
    tris[idx].Z[2] = (pt[2]-10)*0.05;
    tris[idx].normals[2].X = normals[3*ptIds[2]+0];
    tris[idx].normals[2].Y = normals[3*ptIds[2]+1];
    tris[idx].normals[2].Z = normals[3*ptIds[2]+2];

    // 1->2 interpolate between light blue, dark blue
    // 2->2.5 interpolate between dark blue, cyan
    // 2.5->3 interpolate between cyan, green
    // 3->3.5 interpolate between green, yellow
    // 3.5->4 interpolate between yellow, orange
    // 4->5 interpolate between orange, brick
    // 5->6 interpolate between brick, salmon
    double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
    double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
    unsigned char RGB[8][3] = { { 71, 71, 219 }, 
      { 0, 0, 91 },
      { 0, 255, 255 },
      { 0, 128, 0 },
      { 255, 255, 0 },
      { 255, 96, 0 },
      { 107, 0, 0 },
      { 224, 76, 76 } 
    };
    for (int j = 0 ; j < 3 ; j++)
    {
      float val = color_ptr[ptIds[j]];
      int r;
      for (r = 0 ; r < 7 ; r++)
      {
        if (mins[r] <= val && val < maxs[r])
          break;
      }
      if (r == 7)
      {
        cerr << "Could not interpolate color for " << val << endl;
        exit(EXIT_FAILURE);
      }
      double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
      tris[idx].colors[j].R = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
      tris[idx].colors[j].G = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
      tris[idx].colors[j].B = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
    }
  }

  return tris;
}
#endif

void render_flat_bot(Triangle tri, Screen screen, int min, int med, int max) {

  int L;
  int R;

  L = tri.X[min] < tri.X[med] ? min : med;
  R = tri.X[min] < tri.X[med] ? med : min;

  Line *lLine = new Line(tri.X[L], tri.Y[L], tri.X[max], tri.Y[max]);
  Line *rLine = new Line(tri.X[R], tri.Y[R], tri.X[max], tri.Y[max]);

  double rowMin = ceil441(tri.Y[min]);
  double rowMax = floor441(tri.Y[max]);

  for (int r = rowMin; r <= rowMax; r++) {
    double lEnd = lLine->getX(r);
    double rEnd = rLine->getX(r);

    // Get left and right color
    Color lEnd_color = lerp_color(r, tri.Y[L], tri.Y[max], tri.colors[L], tri.colors[max]);
    Color rEnd_color = lerp_color(r, tri.Y[R], tri.Y[max], tri.colors[R], tri.colors[max]);

    // Get left and right normal 
    Normal lEnd_normal = lerp_normal(r, tri.Y[L], tri.Y[max], tri.normals[L], tri.normals[max]);
    Normal rEnd_normal = lerp_normal(r, tri.Y[R], tri.Y[max], tri.normals[R], tri.normals[max]);

    // Get left and right Z
    double lZ = lerp(r, tri.Y[L], tri.Y[max], tri.Z[L], tri.Z[max]);
    double rZ = lerp(r, tri.Y[R], tri.Y[max], tri.Z[R], tri.Z[max]);

    for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
      Color color   = lerp_color(c, lEnd, rEnd, lEnd_color, rEnd_color);
      Normal normal = lerp_normal(c, lEnd, rEnd, lEnd_normal, rEnd_normal);
      double z      = lerp(c, lEnd, rEnd, lZ, rZ);

      screen.ImageColor(r, c, z, color, normal);
    }
  }

  delete lLine;
  delete rLine;
}

void render_flat_top(Triangle tri, Screen screen, int min, int med, int max) {

  int L;
  int R;


  L = tri.X[med] < tri.X[max] ? med : max;
  R = tri.X[med] < tri.X[max] ? max : med;

  //printf("L: %d, R: %d, B: %d\n", L, R, min);

  Line *lLine = new Line(tri.X[min], tri.Y[min], tri.X[L], tri.Y[L]);
  Line *rLine = new Line(tri.X[min], tri.Y[min], tri.X[R], tri.Y[R]);

  double rowMin = ceil441(tri.Y[min]);
  double rowMax = floor441(tri.Y[max]);

  //printf("rowMin: %f, rowMax: %f\n", rowMin, rowMax);

  for (int r = rowMin; r <= rowMax; r++) {
    double lEnd = lLine->getX(r);
    double rEnd = rLine->getX(r);

    //printf("\nr: %d\n", r);
    //printf("  lEnd: %f, rEnd: %f\n", lEnd, rEnd);
    

    // Get color of left and right end
    Color lEnd_color = lerp_color(r, tri.Y[L], tri.Y[min], tri.colors[L], tri.colors[min]);
    Color rEnd_color = lerp_color(r, tri.Y[R], tri.Y[min], tri.colors[R], tri.colors[min]);

    //printf("  lEnd_color: %f, %f, %f rEnd_color: %f, %f, %f\n", 
    //    lEnd_color.R, lEnd_color.G, lEnd_color.B, rEnd_color.R, rEnd_color.G, rEnd_color.B);

    // Get normal of left and right end
    Normal lEnd_normal = lerp_normal(r, tri.Y[L], tri.Y[min], tri.normals[L], tri.normals[min]);
    Normal rEnd_normal = lerp_normal(r, tri.Y[R], tri.Y[min], tri.normals[R], tri.normals[min]);

    //printf("  lEnd_normal: %f, %f, %f rEnd_normal: %f, %f, %f\n", 
    //    lEnd_normal.X, lEnd_normal.Y, lEnd_normal.Z, rEnd_normal.X, rEnd_normal.Y, rEnd_normal.Z);

    // Get left and right Z
    double lZ = lerp(r, tri.Y[L], tri.Y[min], tri.Z[L], tri.Z[min]);
    double rZ = lerp(r, tri.Y[R], tri.Y[min], tri.Z[R], tri.Z[min]);

    for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
      Color color   = lerp_color(c, lEnd, rEnd, lEnd_color, rEnd_color);
      Normal normal = lerp_normal(c, lEnd, rEnd, lEnd_normal, rEnd_normal);
      double z      = lerp(c, lEnd, rEnd, lZ, rZ);

      //printf("\n  c: %d\n", c);
      //printf("    color: %f %f %f\n", color.R, color.G, color.B);
      //printf("    normal: %f %f %f\n", normal.X, normal.Y, normal.Z);
      //printf("    z: %f\n", z);

      //exit(0);
      screen.ImageColor(r, c, z, color, normal);
    }
  }

  delete lLine;
  delete rLine;
}

Triangle split_bot(Triangle tri, int min, int med, int max) {

  // Calculate the line from min to max
  Line *mLine = new Line(tri.X[min], tri.Y[min], tri.X[max], tri.Y[max]);
  Triangle r_tri = tri;

  // Get the x and y values of the new point
  double r_x = mLine->getX(tri.Y[med]);
  double r_y = tri.Y[med];

  // Set the new color, Z, and normals
  r_tri.normals[min] = lerp_normal(r_y, tri.Y[min], tri.Y[max], tri.normals[min], tri.normals[max]);
  r_tri.colors[min]  = lerp_color(r_y, tri.Y[min], tri.Y[max], tri.colors[min], tri.colors[max]);
  r_tri.Z[min]       = lerp(r_y, tri.Y[min], tri.Y[max], tri.Z[min], tri.Z[max]);

  // Set the new X and Y
  r_tri.X[min] = r_x;
  r_tri.Y[min] = r_y;

  delete mLine;
  return r_tri;
}

Triangle split_top(Triangle tri, int min, int med, int max) {

  // Calculate the line from min to max
  Line *mLine = new Line(tri.X[min], tri.Y[min], tri.X[max], tri.Y[max]);

  // Get the x and y values of the new point
  double r_x = mLine->getX(tri.Y[med]);
  double r_y = tri.Y[med];

  Triangle r_tri = tri;

  // Set the new color, Z, and normals
  r_tri.normals[max] = lerp_normal(r_y, tri.Y[min], tri.Y[max], tri.normals[min], tri.normals[max]);
  r_tri.colors[max]  = lerp_color(r_y, tri.Y[min], tri.Y[max], tri.colors[min], tri.colors[max]);
  r_tri.Z[max]       = lerp(r_y, tri.Y[min], tri.Y[max], tri.Z[min], tri.Z[max]);

  r_tri.X[max] = r_x;
  r_tri.Y[max] = r_y;

  delete mLine;
  return r_tri;
}

void render_triangle(Triangle tri, Screen screen) {

  // if horizontally or vertically colinear, skip
  if (tri.Y[0] == tri.Y[1] && tri.Y[1] == tri.Y[2]) return;
  if (tri.X[0] == tri.X[1] && tri.X[1] == tri.X[2]) return;

  // if any two poitns are equal, skip
  if (tri.Y[0] == tri.Y[1] && tri.X[0] == tri.X[1]) return;
  if (tri.Y[1] == tri.Y[2] && tri.X[1] == tri.X[2]) return;
  if (tri.Y[2] == tri.Y[0] && tri.X[2] == tri.X[0]) return;

  // sort the points of the triangle based on Y values
  int min, med, max;

  if( tri.Y[0] > tri.Y[1] ){
    if( tri.Y[0] > tri.Y[2] ){
      max = 0; 
      if( tri.Y[1] > tri.Y[2] ){
        med = 1; 
        min = 2; 
      }else{
        med = 2; 
        min = 1; 
      }
    }else{
      med = 0; 
      if( tri.Y[1] > tri.Y[2] ){
        max = 1; 
        min = 2; 
      }else{
        max = 2; 
        min = 1; 
      }
    }
  }else{
    if( tri.Y[1] > tri.Y[2] ){
      max = 1; 
      if( tri.Y[0] > tri.Y[2] ){
        med = 0; 
        min = 2; 
      }else{
        med = 2;
        min = 0;
      }
    }else{
      med = 1; 
      max = 2;
      min = 0;
    }
  }

  // Handle flat bottom and flat top triangles
  if (tri.Y[min] == tri.Y[med]) {
    render_flat_bot(tri, screen, min, med, max);
    return;
  }
  if (tri.Y[med] == tri.Y[max]) {
    render_flat_top(tri, screen, min, med, max);
    return;
  }


  // Split non-flat triangles
  Triangle tri_bot = split_bot(tri, min, med, max);
  Triangle tri_top = split_top(tri, min, med, max);

  render_flat_bot(tri_bot, screen, min, med, max);
  render_flat_top(tri_top, screen, min, med, max); 

  return;
}

int main()
{
  vtkImageData *image = NewImage(1000, 1000);
  unsigned char *buffer = 
    (unsigned char *) image->GetScalarPointer(0,0,0);
  int npixels = 1000*1000;
  for (int i = 0 ; i < npixels*3 ; i++)
    buffer[i] = 0;

  std::vector<Triangle> triangles = GetTriangles();

  Screen screen;
  screen.buffer = buffer;
  screen.width = 1000;
  screen.height = 1000;
  screen.Zgrid = (double *) malloc(sizeof(double) * screen.width * screen.height);
  for (int i = 0; i < 1000*1000; i++) screen.Zgrid[i] = -1;

  for (int i = 0; i < triangles.size(); i++) {
    //printf("Processing Triangle: %d\n", i);

    render_triangle(triangles[i], screen);
  }

  WriteImage(image, "allTriangles");
}
