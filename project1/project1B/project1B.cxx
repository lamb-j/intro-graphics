#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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
    double         X[3];
    double         Y[3];
    unsigned char color[3];
    int T, L, R;

    Line *lLine;
    Line *rLine;

    // would some methods for transforming the triangle in place be helpful?
    int getT();
    int getL();
    int getR();
    void setPoints();
    void createLines();
    void freeLines();
};

int Triangle::getT() {
  if (Y[0] > Y[1] && Y[0] > Y[2]) return 0;
  if (Y[1] > Y[0] && Y[1] > Y[2]) return 1;
  if (Y[2] > Y[1] && Y[2] > Y[0]) return 2;

  return -1;
}

int Triangle::getL() {
  if (T == -1) return -1;
 
  if (T != 0 && X[0] <= X[1] && X[0] <= X[2]) return 0;
  if (T != 1 && X[1] <= X[0] && X[1] <= X[2]) return 1;
  if (T != 2 && X[2] <= X[0] && X[2] <= X[1]) return 2;

  return -1;
}

int Triangle::getR() {
  if (T == -1 || L == -1) return -1;

  for (int i = 0; i < 3; i++)
    if (T != i && L != i) return i;

  return -1;
}

void Triangle::createLines() {
  setPoints();

  lLine = new Line(X[T], Y[T], X[L], Y[L]);
  rLine = new Line(X[T], Y[T], X[R], Y[R]);
}

void Triangle::freeLines() {
  delete lLine;
  delete rLine;
}

void Triangle::setPoints() {
  T = getT();
  L = getL();
  R = getR();
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

std::vector<Triangle> GetTriangles(void)
{
  std::vector<Triangle> rv(100);

  unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
    {76,153,0}, {255, 204, 204}, {204, 204, 0}};
  for (int i = 0 ; i < 100 ; i++)
  {
    int idxI = i%10;
    int posI = idxI*100;
    int idxJ = i/10;
    int posJ = idxJ*100;
    int firstPt = (i%3);
    rv[i].X[firstPt] = posI;
    if (i == 50)
      rv[i].X[firstPt] = -10;
    rv[i].Y[firstPt] = posJ;
    rv[i].X[(firstPt+1)%3] = posI+99;
    rv[i].Y[(firstPt+1)%3] = posJ;
    rv[i].X[(firstPt+2)%3] = posI+i;
    rv[i].Y[(firstPt+2)%3] = posJ+10*(idxJ+1);
    if (i == 95)
      rv[i].Y[(firstPt+2)%3] = 1050;
    rv[i].color[0] = colors[i%6][0];
    rv[i].color[1] = colors[i%6][1];
    rv[i].color[2] = colors[i%6][2];
  }

  return rv;
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

  // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
  for (int i = 0; i < triangles.size(); i++) {
    //printf("Processing: Triangle %d\n", i);

    triangles[i].createLines();

    int L = triangles[i].L;
    int T = triangles[i].T;
    double rowMin = ceil441(triangles[i].Y[L]);
    double rowMax = floor441(triangles[i].Y[T]);

    // Scanline for this triangle
    for (int r = rowMin; r <= rowMax; r++) {
      double lEnd = triangles[i].lLine->getX(r);
      double rEnd = triangles[i].rLine->getX(r);

      for (int c = ceil441(lEnd); c <= floor441(rEnd); c++) {
        screen.ImageColor(r, c, triangles[i].color);
      }

    }

    triangles[i].freeLines();
  }

  WriteImage(image, "allTriangles");
}
