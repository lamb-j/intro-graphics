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
#include <cmath>

using std::cerr;
using std::endl;

class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

  void
Matrix::Print(ostream &o)
{
  for (int i = 0 ; i < 4 ; i++)
  {
    char str[256];
    sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
    o << str;
  }
}

  Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
  Matrix rv;
  for (int i = 0 ; i < 4 ; i++)
    for (int j = 0 ; j < 4 ; j++)
    {
      rv.A[i][j] = 0;
      for (int k = 0 ; k < 4 ; k++)
        rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
    }

  return rv;
}

  void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
  ptOut[0] = ptIn[0]*A[0][0]
    + ptIn[1]*A[1][0]
    + ptIn[2]*A[2][0]
    + ptIn[3]*A[3][0];
  ptOut[1] = ptIn[0]*A[0][1]
    + ptIn[1]*A[1][1]
    + ptIn[2]*A[2][1]
    + ptIn[3]*A[3][1];
  ptOut[2] = ptIn[0]*A[0][2]
    + ptIn[1]*A[1][2]
    + ptIn[2]*A[2][2]
    + ptIn[3]*A[3][2];
  ptOut[3] = ptIn[0]*A[0][3]
    + ptIn[1]*A[1][3]
    + ptIn[2]*A[2][3]
    + ptIn[3]*A[3][3];
}

struct LightingParameters
{
  LightingParameters(void)
  {
    lightDir[0] = -0.6;
    lightDir[1] = 0;
    lightDir[2] = -0.8;
    Ka = 0.3;
    Kd = 0.7;
    Ks = 5.3;
    alpha = 7.5;
  };

  double lightDir[3];  // The direction of the light source
  double Ka;           // The coefficient for ambient lighting.
  double Kd;           // The coefficient for diffuse lighting.
  double Ks;           // The coefficient for specular lighting.
  double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;
class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    //Matrix          ViewTransform(void) {;};
    //Matrix          CameraTransform(void) {;};
    //Matrix          DeviceTransform(void) {;};
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
  int nNonRamp = nFrames-2*ramp;
  double height = 1./(nNonRamp + 4*ramp/M_PI);
  if (curFrame < ramp)
  {
    double factor = 2*height*ramp/M_PI;
    double eval = cos(M_PI/2*((double)curFrame)/ramp);
    return (1.-eval)*factor;
  }
  else if (curFrame > nFrames-ramp)
  {
    int amount_left = nFrames-curFrame;
    double factor = 2*height*ramp/M_PI;
    double eval =cos(M_PI/2*((double)amount_left/ramp));
    return 1. - (1-eval)*factor;
  }
  double amount_in_quad = ((double)curFrame-ramp);
  double quad_part = amount_in_quad*height;
  double curve_part = height*(2*ramp)/M_PI;
  return quad_part+curve_part;
}

  Camera
GetCamera(int frame, int nframes)
{
  double t = SineParameterize(frame, nframes, nframes/10);
  Camera c;
  c.near = 5;
  c.far = 200;
  c.angle = M_PI/6;
  c.position[0] = 40*sin(2*M_PI*t);
  c.position[1] = 40*cos(2*M_PI*t);
  c.position[2] = 40;
  c.focus[0] = 0;
  c.focus[1] = 0;
  c.focus[2] = 0;
  c.up[0] = 0;
  c.up[1] = 1;
  c.up[2] = 0;
  return c;
}

double deg_2_rad(double deg) {
  return (deg * M_PI) / 180.0;
}


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

// Technically not all normals, but I didn't want to name the class Vector...
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

    Normal operator+(const Normal& b) {
      Normal normal;
      normal.X = this->X + b.X;
      normal.Y = this->Y + b.Y;
      normal.Z = this->Z + b.Z;

      return normal;
    }

    Normal operator-(const Normal& b) {
      Normal normal;
      normal.X = this->X - b.X;
      normal.Y = this->Y - b.Y;
      normal.Z = this->Z - b.Z;

      return normal;
    }

    Normal operator*(double d) {
      Normal n;
      n.X = this->X * d;
      n.Y = this->Y * d;
      n.Z = this->Z * d;

      return n;
    }
};

Normal cross_product(Normal a, Normal b) {
  Normal c;

  c.X = a.Y * b.Z - a.Z * b.Y;
  c.Y = a.Z * b.X - a.X * b.Z;
  c.Z = a.X * b.Y - a.Y * b.X;

  return c;
}

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
    double  shade[3];

    void print_triangle();
};

void Triangle::print_triangle() {
  for (int i = 0; i < 3; i++) {
    printf("Tri %d: X: %f, Y: %f, Z: %f\n R: %f, G: %f, B: %f\n, Nx: %f, Ny: %f, Nz: %f\n", 
        i, X[i], Y[i], Z[i], colors[i].R, colors[i].G, colors[i].B, normals[i].X, normals[i].Y, normals[i].Z);
  }
}
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
  /*
     vtkDataSetWriter *writer = vtkDataSetWriter::New();
     writer->SetInput(pd);
     writer->SetFileName("hrc.vtk");
     writer->Write();
   */

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
    tris[idx].X[0] = pt[0];
    tris[idx].Y[0] = pt[1];
    tris[idx].Z[0] = pt[2];
    tris[idx].normals[0].X = normals[3*ptIds[0]+0];
    tris[idx].normals[0].Y = normals[3*ptIds[0]+1];
    tris[idx].normals[0].Z = normals[3*ptIds[0]+2];
    pt = pts->GetPoint(ptIds[1]);
    tris[idx].X[1] = pt[0];
    tris[idx].Y[1] = pt[1];
    tris[idx].Z[1] = pt[2];
    tris[idx].normals[1].X = normals[3*ptIds[1]+0];
    tris[idx].normals[1].Y = normals[3*ptIds[1]+1];
    tris[idx].normals[1].Z = normals[3*ptIds[1]+2];
    pt = pts->GetPoint(ptIds[2]);
    tris[idx].X[2] = pt[0];
    tris[idx].Y[2] = pt[1];
    tris[idx].Z[2] = pt[2];
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

double CalculatePhongShading(LightingParameters &, Normal normal, Normal viewDir)
{
  Normal light;
  light.X = lp.lightDir[0];
  light.Y = lp.lightDir[1];
  light.Z = lp.lightDir[2];

  double dp = dot_product(light.normalize(), normal.normalize());

  // Ambient switch
  double ambient = 1;

  // Calulate two-sided for diffuse
  double diffuse = std::abs(dp); 

  // Calculate one-sided for specular
  Normal R; 
  R.X = 2 * dp * normal.X - light.X;
  R.Y = 2 * dp * normal.Y - light.Y;
  R.Z = 2 * dp * normal.Z - light.Z;

  double specular = pow(std::max(0.0, dot_product(viewDir.normalize(), R.normalize())), lp.alpha);

  return lp.Ka*ambient + lp.Kd*diffuse + lp.Ks*specular;
}

class Screen
{
  public:
    unsigned char *buffer;
    double *Zgrid;
    int width, height;
    vtkImageData *image;

    void ImageColor(int r, int c, double z, Color color, double shade);
    void Allocate(int width, int height);
    void Initalize();
    void Free();
};

void Screen::Free() {
  delete buffer;
  delete Zgrid;
}

// Make sure you allocate first
void Screen::Initalize() {
  for (int i = 0; i < width*height; i++) {
    Zgrid[i] = -2;
  }
  for (int i = 0; i < 3*width*height; i++) {
    buffer[i] = 0;
  }
}

void Screen::Allocate(int width, int height) {
  this->width = width;
  this->height = height;

  this->image = NewImage(width, height); 
  this->buffer = 
    (unsigned char *) image->GetScalarPointer(0,0,0);

  Zgrid = (double *) malloc(sizeof(double) * width * height);
}

void Screen::ImageColor(int r, int c, double z, Color color, double shade) {
  if (r >= height || r < 0 || c >= width || c < 0) return;

  if (z < Zgrid[r*width + c]) return;

  Zgrid[r*width + c] = z;

  Normal V;
  V.X = 0;
  V.Y = 0;
  V.Z = -1;

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

    // Get left and right Z
    double lZ = lerp(r, tri.Y[L], tri.Y[max], tri.Z[L], tri.Z[max]);
    double rZ = lerp(r, tri.Y[R], tri.Y[max], tri.Z[R], tri.Z[max]);

    // Get left and right shade
    double lShade = lerp(r, tri.Y[L], tri.Y[max], tri.shade[L], tri.shade[max]);
    double rShade = lerp(r, tri.Y[R], tri.Y[max], tri.shade[R], tri.shade[max]);

    for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
      Color color   = lerp_color(c, lEnd, rEnd, lEnd_color, rEnd_color);
      double z      = lerp(c, lEnd, rEnd, lZ, rZ);
      double shade  = lerp(c, lEnd, rEnd, lShade, rShade); 

      screen.ImageColor(r, c, z, color, shade);
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

  Line *lLine = new Line(tri.X[min], tri.Y[min], tri.X[L], tri.Y[L]);
  Line *rLine = new Line(tri.X[min], tri.Y[min], tri.X[R], tri.Y[R]);

  double rowMin = ceil441(tri.Y[min]);
  double rowMax = floor441(tri.Y[max]);

  for (int r = rowMin; r <= rowMax; r++) {
    double lEnd = lLine->getX(r);
    double rEnd = rLine->getX(r);

    // Get color of left and right end
    Color lEnd_color = lerp_color(r, tri.Y[L], tri.Y[min], tri.colors[L], tri.colors[min]);
    Color rEnd_color = lerp_color(r, tri.Y[R], tri.Y[min], tri.colors[R], tri.colors[min]);

    // Get left and right Z
    double lZ = lerp(r, tri.Y[L], tri.Y[min], tri.Z[L], tri.Z[min]);
    double rZ = lerp(r, tri.Y[R], tri.Y[min], tri.Z[R], tri.Z[min]);

    // Get left and right shade
    double lShade = lerp(r, tri.Y[L], tri.Y[min], tri.shade[L], tri.shade[min]);
    double rShade = lerp(r, tri.Y[R], tri.Y[min], tri.shade[R], tri.shade[min]);

    for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
      Color color   = lerp_color(c, lEnd, rEnd, lEnd_color, rEnd_color);
      double z      = lerp(c, lEnd, rEnd, lZ, rZ);
      double shade  = lerp(c, lEnd, rEnd, lShade, rShade);

      screen.ImageColor(r, c, z, color, shade);
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

  // Set the new color, Z, normals, and shade
  r_tri.colors[min]  = lerp_color(r_y, tri.Y[min], tri.Y[max], tri.colors[min], tri.colors[max]);
  r_tri.Z[min]       = lerp(r_y, tri.Y[min], tri.Y[max], tri.Z[min], tri.Z[max]);
  r_tri.shade[min]   = lerp(r_y, tri.Y[min], tri.Y[max], tri.shade[min], tri.shade[max]);

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

  // Set the new color, Z, normals, and shade
  r_tri.colors[max]  = lerp_color(r_y, tri.Y[min], tri.Y[max], tri.colors[min], tri.colors[max]);
  r_tri.Z[max]       = lerp(r_y, tri.Y[min], tri.Y[max], tri.Z[min], tri.Z[max]);
  r_tri.shade[max]   = lerp(r_y, tri.Y[min], tri.Y[max], tri.shade[min], tri.shade[max]);

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

Triangle CalculateShade(Triangle triangle, Camera c) {

  for (int i = 0; i < 3; i++) {
    Normal viewDir;
    viewDir.X = triangle.X[i] - c.position[0];
    viewDir.Y = triangle.Y[i] - c.position[1];
    viewDir.Z = triangle.Z[i] - c.position[2];

    viewDir = viewDir.normalize();

    triangle.shade[i] = CalculatePhongShading(lp, triangle.normals[i].normalize(), viewDir.normalize());
  }

  return triangle;
}

Triangle TransformTriangle(Triangle t, Matrix Compose) {

  // Transform each point
  for (int i = 0; i < 3; i++) {
    double pIn[4];
    double pOut[4];

    pIn[0] = t.X[i];
    pIn[1] = t.Y[i];
    pIn[2] = t.Z[i];
    pIn[3] = 1;

    Compose.TransformPoint(pIn, pOut);

    // divide all by w
    for (int j = 0; j < 3; j++) pOut[j] /= pOut[3];

    t.X[i] = pOut[0];
    t.Y[i] = pOut[1];
    t.Z[i] = pOut[2];
  }
  return t;
}

int main()
{
  std::vector<Triangle> triangles = GetTriangles();

  Screen screen;
  screen.Allocate(1000, 1000);

  for (int i = 0; i < 4; i++) {
    printf("Processing Camera %d, %d\n", i, 1000);

    screen.Initalize();

    // Set up camera
    Camera c = GetCamera(i*250, 1000);

    // Calcualte shading at each triangle
    for (int j = 0; j < triangles.size(); j++) {
      triangles[j] = CalculateShade(triangles[j], c);
    }

    // Calculate camera frame
    Normal O;
    O.X = c.position[0];
    O.Y = c.position[1];
    O.Z = c.position[2];

    Normal focus;
    focus.X = c.focus[0];
    focus.Y = c.focus[1];
    focus.Z = c.focus[2];

    Normal Up;
    Up.X = c.up[0];
    Up.Y = c.up[1];
    Up.Z = c.up[2];

    Normal v3 = O - focus;
    Normal v1 = cross_product(Up, v3);
    Normal v2 = cross_product(v3, v1);

    v1 = v1.normalize();
    v2 = v2.normalize();
    v3 = v3.normalize();

    Normal t;
    t.X = 0 - O.X;
    t.Y = 0 - O.Y;
    t.Z = 0 - O.Z;

    Matrix CameraTransform;
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        CameraTransform.A[j][k] = 0;
      }
    }

    CameraTransform.A[0][0] = v1.X;
    CameraTransform.A[0][1] = v2.X;
    CameraTransform.A[0][2] = v3.X;
    CameraTransform.A[0][3] = 0; 

    CameraTransform.A[1][0] = v1.Y;
    CameraTransform.A[1][1] = v2.Y;
    CameraTransform.A[1][2] = v3.Y;
    CameraTransform.A[1][3] = 0;

    CameraTransform.A[2][0] = v1.Z;
    CameraTransform.A[2][1] = v2.Z;
    CameraTransform.A[2][2] = v3.Z;
    CameraTransform.A[2][3] = 0;

    CameraTransform.A[3][0] = dot_product(v1, t);
    CameraTransform.A[3][1] = dot_product(v2, t);
    CameraTransform.A[3][2] = dot_product(v3, t);
    CameraTransform.A[3][3] = 1;

    Matrix ViewTransform;
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        ViewTransform.A[j][k] = 0;
      }
    }

    ViewTransform.A[0][0] = 1 / (tan(c.angle / 2));
    ViewTransform.A[1][1] = 1 / (tan(c.angle / 2));
    ViewTransform.A[2][2] = (c.far + c.near) / (c.far - c.near);
    ViewTransform.A[2][3] = -1;
    ViewTransform.A[3][2] = (2 * c.far * c.near) / (c.far - c.near);

    Matrix DeviceTransform;
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        DeviceTransform.A[j][k] = 0;
      }
    }

    double W = screen.width / 2;
    DeviceTransform.A[0][0] = W;
    DeviceTransform.A[1][1] = W;
    DeviceTransform.A[2][2] = 1;

    DeviceTransform.A[3][0] = W;
    DeviceTransform.A[3][1] = W; 
    DeviceTransform.A[3][2] = 0;
    DeviceTransform.A[3][3] = 1;

    Matrix Compose;
    Compose = Matrix::ComposeMatrices(Matrix::ComposeMatrices(CameraTransform, ViewTransform), DeviceTransform);

    // Render triangles 
    for (int j = 0; j < triangles.size(); j++) {

      // Transform the triangle
      Triangle t = TransformTriangle(triangles[j], Compose);
      render_triangle(t, screen);
    }

    // save image
    char buffer[50];
    sprintf(buffer, "frame%d", i*250);
    WriteImage(screen.image, buffer);

  }
  screen.Free();
}



