// Jacob Lambert
// CIS 541
// Project 1A
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1024, 1350); 
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);

   for (int i = 0; i < 27; i++) {
     unsigned char R, G, B;

     // Set the color for this stripe
     switch(i % 3) {
       case 0: B = 0; break;
       case 1: B = 128; break;
       case 2: B = 255; break;
     }

     switch((i / 3) % 3) {
       case 0: G = 0; break;
       case 1: G = 128; break;
       case 2: G = 255; break;
     }

     switch(i / 9) {
       case 0: R = 0; break;
       case 1: R = 128; break;
       case 2: R = 255; break;
     }
     
     // Determine where the stripe should start
     int offset = i * 50 * (3 * 1024);

     // Paint each pixel in the stripe
     for (int j = 0; j < 1024*50; j++) {
       buffer[offset + 3*j] = R;
       buffer[offset + 3*j + 1] = G;
       buffer[offset + 3*j + 2] = B;
     }
  }
  
   WriteImage(image, "allColors");
}
