#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
// includes from get_triangle.cxx
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

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
    double         X[3];
    double         Y[3];
    unsigned char color[3];

    void print_triangle();
};

void Triangle::print_triangle() {
  printf("tri.X[0] = %f\n", X[0]);
  printf("tri.X[1] = %f\n", X[1]);
  printf("tri.X[2] = %f\n\n", X[2]);

  printf("tri.Y[0] = %f\n", Y[0]);
  printf("tri.Y[1] = %f\n", Y[1]);
  printf("tri.Y[2] = %f\n\n", Y[2]);

  printf("tri.color[0] = %d\n", color[0]);
  printf("tri.color[1] = %d\n", color[1]);
  printf("tri.color[2] = %d\n", color[2]);
}

class Screen
{
  public:
    unsigned char   *buffer;
    int width, height;


    // would some methods for accessing and setting pixels be helpful?
    void ImageColor(int r, int c, unsigned char color[3]);
};

void Screen::ImageColor(int r, int c, unsigned char color[3]) {
  if (r >= height || r < 0 || c >= width || c < 0) return;

  int offset = 3*(r*width + c);
  buffer[offset]     = color[0];
  buffer[offset + 1] = color[1];
  buffer[offset + 2] = color[2];
}
  std::vector<Triangle>
GetTriangles(void)
{
  vtkPolyDataReader *rdr = vtkPolyDataReader::New();
  rdr->SetFileName("proj1c_geometry.vtk");
  cerr << "Reading" << endl;
  rdr->Update();
  if (rdr->GetOutput()->GetNumberOfCells() == 0)
  {
    cerr << "Unable to open file!!" << endl;
    exit(EXIT_FAILURE);
  }
  vtkPolyData *pd = rdr->GetOutput();
  int numTris = pd->GetNumberOfCells();
  vtkPoints *pts = pd->GetPoints();
  vtkCellArray *cells = pd->GetPolys();
  vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
  float *color_ptr = colors->GetPointer(0);
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
    tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
    tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
    tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
    tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
    tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
    tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
    tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
    tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
    tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
  }
  cerr << "Done reading" << endl;

  return tris;
}

#if 0
std::vector<Triangle> GetTriangles(void)
{
  std::vector<Triangle> rv(1);

  rv[0].X[0] = 100;
  rv[0].X[1] = 100;
  rv[0].X[2] = -400;

  rv[0].Y[0] = 100;
  rv[0].Y[1] = 800;
  rv[0].Y[2] = 400;

  rv[0].color[0] = 255;
  rv[0].color[1] = 0;
  rv[0].color[2] = 0;

  return rv;
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

    for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
      screen.ImageColor(r, c, tri.color);
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

    for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
      screen.ImageColor(r, c, tri.color);
    }
  }

  delete lLine;
  delete rLine;
}

Triangle split_bot(Triangle tri, int min, int med, int max) {

  // Calculate the line from min to max
  Line *mLine = new Line(tri.X[min], tri.Y[min], tri.X[max], tri.Y[max]);

  // Get the x and y values of the new point
  double r_x = mLine->getX(tri.Y[med]);
  double r_y = tri.Y[med];

  tri.X[min] = r_x;
  tri.Y[min] = r_y;

  delete mLine;
  return tri;
}

Triangle split_top(Triangle tri, int min, int med, int max) {

  // Calculate the line from min to max
  Line *mLine = new Line(tri.X[min], tri.Y[min], tri.X[max], tri.Y[max]);

  // Get the x and y values of the new point
  double r_x = mLine->getX(tri.Y[med]);
  double r_y = tri.Y[med];

  tri.X[max] = r_x;
  tri.Y[max] = r_y;

  delete mLine;
  return tri;
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
  vtkImageData *image = NewImage(1786, 1344);
  unsigned char *buffer = 
    (unsigned char *) image->GetScalarPointer(0,0,0);
  int npixels = 1786*1344;
  for (int i = 0 ; i < npixels*3 ; i++)
    buffer[i] = 0;

  std::vector<Triangle> triangles = GetTriangles();

  Screen screen;
  screen.buffer = buffer;
  screen.width = 1786;
  screen.height = 1344;

  for (int i = 0; i < triangles.size(); i++) {
    render_triangle(triangles[i], screen);
  }

  WriteImage(image, "allTriangles");
}
