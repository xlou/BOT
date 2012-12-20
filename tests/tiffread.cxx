//////////////////////////////////////////////////////////////////////////////////////////////
/*
Author: S.astien Tosi
Institution: IRB Barcelona
Date: 18/12/2012

Description:
Very simple example showing the usage of libtiff to read a (multi-page) tiff image from file

Limitations: 
- The tiff file is not striped nor tiled
- The images within a stack all have same size
- Supported formats: 8-bit greyscale or 32-bit RGB
*/
///////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
// <stdio.h> must be put before others, see
// http://stackoverflow.com/questions/6678092/how-to-resolve-error-c2664-vswprintf-c-l-error-in-visual-studio-2005
#include "tiffio.h"
#include <iostream>
#include "../include/TIFFReaderWriter.hxx"

using namespace std;

int main(int argc, char* argv[])
{
	tsize_t row;	// Image height
	tsize_t col;	// Image width
	tsize_t bytes;	// Number of bytes per scanline
	short bps;		// Number of bit per pixel
	TIFF* tifPtr;	// Pointer to TIFF structure
	unsigned char* ImageBuffer=NULL;	// Pointer to image buffer
	uint16* ImageBuffer16=NULL;
	uint32* RGBImageBuffer = NULL;		// Pointer to RGB image buffer
	int i,j,Planes=0;

	// Open image from file
	tifPtr = TIFFOpen(argv[1], "r");
	
	// Access data
	if (tifPtr!=NULL)
	{
        do 
		{
            Planes++;
			TIFFGetField(tifPtr, TIFFTAG_IMAGELENGTH, &row);
			TIFFGetField(tifPtr, TIFFTAG_IMAGEWIDTH, &col);
			TIFFGetField(tifPtr,TIFFTAG_BITSPERSAMPLE, &bps);
			bytes = TIFFScanlineSize(tifPtr); 
			printf("\n\nImage %u size: %u x %u (%u), bps=%u\n\n", Planes, row, col, bytes, bps);

			// Only handle 8 bit image (greyscale or RGB)
			if(bps==8)
			{
				if(col==bytes) {
					if(Planes==1)ImageBuffer = static_cast<unsigned char* >( _TIFFmalloc(row*col) );
					if(ImageBuffer != NULL)for (i = 0; i < row; i++)TIFFReadScanline(tifPtr, ImageBuffer+i*bytes, i, 0);
				}
				else {
					if(Planes==1)RGBImageBuffer = (uint32*) _TIFFmalloc(row*col* sizeof (uint32));
					if(RGBImageBuffer != NULL)TIFFReadRGBAImage(tifPtr, col, row, RGBImageBuffer , 0);
				}	
			}
			else if (bps == 16) {
/*				std::cerr << (int)row << " " << col << std::endl;
                if (Planes == 1)
					ImageBuffer16 = static_cast<uint16* >( _TIFFmalloc(row * col * sizeof(uint16)) );
                if (ImageBuffer16 != NULL) {
					for (i = 0; i < row; i++)
						TIFFReadScanline(tifPtr, ImageBuffer16+i*row, i, 0);
				}
				else {
					std::cerr << "_TIFFmalloc: fail to allocate memory!!" << std::endl;
					return -1;
				}*/
            }
			else 
				printf("Unsupported format\n\n");
/*
			// Print first row
            if(ImageBuffer != NULL) {
                for (i=0;i<row;i++) {
                    for(j=0;j<col;j++)
                        printf("%u ", ImageBuffer[i*col+j]);
                    printf("\n");
                }
            }

			if(ImageBuffer16 != NULL) {
				for (i=0; i<row; i++) {
	                for(j=0; j<col; j++)
    	                printf("%u ", ImageBuffer16[i*col+j]);
					printf("\n");
				}
			}

			if(RGBImageBuffer!=NULL)
				for(j=0;j<col;j++)
					printf("%x ",RGBImageBuffer[i*0+j]);
*/
		}
		while (TIFFReadDirectory(tifPtr));
		
		if(ImageBuffer!=NULL)_TIFFfree(ImageBuffer);
		if(ImageBuffer16!=NULL)_TIFFfree(ImageBuffer16);
		if(RGBImageBuffer!=NULL)_TIFFfree(RGBImageBuffer);
		TIFFClose(tifPtr);
	}


	// try TIFFReaderWriter
	Matrix2D image = Matrix2D(Matrix2D::difference_type((int)row, (int)col));
	std::cerr << "good here" << std::endl;
	Matrix2D out2d = TIFFReaderWriter::loadTiff(argv[1]);
	std::cout << "\n" << out2d << std::endl;
}
