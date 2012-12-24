/*!
 * @author  Xinghua Lou <xinghua.lou@iwr.uni-heidelberg.de>
 *
 * @section LICENSE
 * 
 * BOT. Copyright (c) 2010 by Xinghua Lou.
 *
 * This software was developed by Xinghua Lou.
 * Enquiries shall be directed to: xinghua.lou@iwr.uni-heidelberg.de.
 *
 * All advertising materials mentioning features or use of this software must
 * display the following acknowledgement: ``This product includes the BOT
 * library developed by Xinghua Lou. Please direct enquiries concerning BOT to 
 * xinghua.lou@iwr.uni-heidelberg.de.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, 
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - All advertising materials mentioning features or use of this software must 
 *   display the following acknowledgement: ``This product includes the BOT
 *   library developed by Xinghua Lou. Please direct enquiries concerning BOT to 
 *   xinghua.lou@iwr.uni-heidelberg.de.
 * - The names of the authors must not be used to endorse or promote products 
 *   derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
 * EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __TIFF_READER_WRITER_HXX
#define __TIFF_READER_WRITER_HXX

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include "tiffio.h"
#include "TypeDefinition.hxx"
#include "SolutionCoder.hxx"
#include "TrainingData.hxx"
#include "InputOutput.hxx"
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include <vigra/matrix.hxx> 


using namespace vigra::linalg;
//using namespace std;

namespace bot 
{

/*! This class allows for reading the images/segmentations from and writing 
 *  the tracking results to a hdf5 file.
 *
 */
class TIFFReaderWriter
{
public:
	static bool loadTiff(const std::string& file, Matrix2D& image)
	{
		tsize_t row_t;    // Image height
	    tsize_t col_t;    // Image width
	    short bps;      // Number of bit per pixel
	    TIFF* tifPtr;   // Pointer to TIFF structure
    	int zInd = 0;

	    // Open image from file
	    tifPtr = TIFFOpen(file.c_str(), "r");

	    // Access data
	    if (tifPtr == NULL) 
			return false;

	    TIFFGetField(tifPtr, TIFFTAG_IMAGELENGTH, &row_t);
	    TIFFGetField(tifPtr, TIFFTAG_IMAGEWIDTH, &col_t);
	    TIFFGetField(tifPtr, TIFFTAG_BITSPERSAMPLE, &bps);
		int row = static_cast<int >(row_t);
		int col = static_cast<int >(col_t);

		// Fixme: now we only read the first layer, 3d tiff stack is not supported!!
        do {
            zInd ++;

			// create output matrix
			image.reshape(Matrix2D::difference_type(row, col));

			switch (bps) {
			case 8:
			{
				// readin line by line
				unsigned char* buf = static_cast<unsigned char* >(_TIFFmalloc(col));
				for (int i = 0; i < row; i++) {
					TIFFReadScanline(tifPtr, buf, i, 0);
					for (int j = 0; j < col; j++)
						image(i, j) = buf[j];
				}
				_TIFFfree(buf);
			}
				break;
			case 16: 
			{
                // readin line by line
                uint16* buf = static_cast<uint16* >(_TIFFmalloc(col * sizeof(uint16)));
                for (int i = 0; i < row; i++) {
                    TIFFReadScanline(tifPtr, buf, i, 0);
                    for (int j = 0; j < col; j++)
                        image(i, j) = buf[j];
				}
				_TIFFfree(buf);
			}
				break;
			default:
				std::cerr << "Unsupproted file format!" << std::endl;
			}
		}
    	while (TIFFReadDirectory(tifPtr) && zInd < 1); // only read the first slice, need to be exteded in the future

		return true;
	}

	static bool loadTiffDir(const std::string& dir,
                   std::vector<Matrix2D >& matrices)
	{
		std::vector<std::string > files;
		files = getFilesInDir(dir, ".tiff");
		if (files.size() == 0)
			 files = getFilesInDir(dir, ".tif");


		for (int i=0; i<files.size(); i++) {
			Matrix2D image;
			loadTiff(files[i], image);
			matrices.push_back(image);
		}

		return true;
	}

    static bool loadTiffDir(const std::string& dir,
					        std::vector<Matrix2D >& matrices, 
							std::vector<std::string >& files)
	{
		files = getFilesInDir(dir, ".tiff");
		if (files.size() == 0)
			files = getFilesInDir(dir, ".tif");

        for (int i=0; i<files.size(); i++) {
			Matrix2D image;
		    loadTiff(files[i], image);
			matrices.push_back(image);
		}
 
		return true;
	}
	
	static bool parseAnnotationLine(const std::string& line, Matrix2D& source, Matrix2D& target)
	{
		int pos = line.find_first_of("->");
		if (pos == std::string::npos)
			return false;

		std::vector<MatrixElem > tokens;

		// parse the part before "->", viz. source
		tokens = VigraSTLInterface::string_to_vector<MatrixElem >(line.substr(0, pos));
		source = VigraSTLInterface::vector_to_matrix<MatrixElem >(tokens);
		
		// parse the part after "->", viz. target
		tokens = VigraSTLInterface::string_to_vector<MatrixElem >(line.substr(pos+2));
		target = VigraSTLInterface::vector_to_matrix<MatrixElem >(tokens);
		
		return true;
	};

    /*! A static function that loads the training data
     *  @param filename The path to the hdf5 file
     *  @param training The training data
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool loadAnnotationDir(const std::string& dir,
		   const std::vector<std::string >& references,
                   TrainingData& training) 
    {
		std::vector<std::string > files = getFilesInDir(dir, ".txt");
		char buf[256];

		for (int i=0; i<files.size(); i++) {
			int time = findRawFileIndex(references, files[i]);
			std::cerr << "Annotation " << files[i] << " -> " << references[time] << std::endl;
			// open file
			std::ifstream ifs(files[i].c_str(), std::ifstream::in);
			
			// create object associations
			LabelAssociations associations;
			
			// variables to store association for a single event
			LabelAssociation association;
			Matrix2D source;
			Matrix2D target;
			std::string name;

			// read line by line
			while (!ifs.eof()) {
				ifs.getline(buf, 256);
				std::string line(buf);
				if (line.find_first_of("[") != std::string::npos && line.find_first_of("]") != std::string::npos) { // indicator of new event 
					// add last association, if necessary
					if (source.shape(0) != 0 && source.shape(0) == target.shape(0)) {
						association.name = name;
						association.source = source;
						association.target = target;
						associations.push_back(association);
					}
				
					// parse event name
					int leftBr = line.find_first_of("[");
					int rightBr = line.find_first_of("]");
					name = line.substr(leftBr+1, rightBr - leftBr - 1);
					
					// reset association object and source/target
					association = LabelAssociation();
					source = Matrix2D();
					target = Matrix2D();
				}
				else if (line.find_first_of("->") != std::string::npos) { // annotated events
					if (source.elementCount() == 0 && target.elementCount() == 0) { // first association for a new event
						// dirstly assign variable source and target
						Matrix2D source_, target_;
						parseAnnotationLine(line, source_, target_);
						source = transpose(source_);
						target = transpose(target_);
					}
					else {
						// cascade to source and target
						Matrix2D source_, target_;
						parseAnnotationLine(line, source_, target_);
						try {
							source = joinVertically(source, transpose(source_));
							target = joinVertically(target, transpose(target_));
						}
						catch (...) {
							std::cerr << "Fail to parse line '" << line << "' for event" << association.name << std::endl;
							return false;
						}
					}
				}
				else { // line with unknown format, skip it
					
				}
			}
			
			// add last association, if necessary
			if (source.shape(0) != 0 && source.shape(0) == target.shape(0)) {
				association.name = name;
				association.source = source;
				association.target = target;
				associations.push_back(association);
			}
			
			training.add(time, associations);
            
            ifs.close();
		}
    };

    /*! A static function that saves the label associations at a given time
     *  @param filename The path to the hdf5 file
     *  @param associations The label associations
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool saveAssociationFile(
        const int time,
        const std::string& file, 
        const LabelAssociations& associations)
    {
        std::cerr << "\tTime=" << time << "; result=" << file << std::endl;
        
		std::ofstream ofs(file.c_str(), ios_base::in | ios_base::out | ios_base::trunc);
		for (int i=0; i<associations.size(); i++) {
			// output event name
			LabelAssociation association = associations[i];
			ofs << "[" << association.name << "]" << std::endl;
			
			// output associations
			Matrix2D source = association.source;
			Matrix2D target = association.target;
			if (source.shape(0) == 0 || target.shape(0) == 0)
				continue ;
			for (int j=0; j<source.shape(0); j++) {
				int k;
				for (k=0; k<source.shape(1); k++) 
					ofs << source(j, k) << " ";
				ofs << "->";
				for (k=0; k<target.shape(1); k++) 
					ofs << " " << target(j, k);
				ofs << std::endl;
			}
			ofs << std::endl;
		}
		
        ofs.close();
		return true;
    };

    /*! A static function that saves all the tracking solutions
     *  @param filename The path to the hdf5 file
     *  @param framepairs The vector of framepairs
     *  @param singlets_vec The vector of Singlets objects
     *  @param multiplets_vec The vector of Multiplets objects
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool saveAssociations(
        const std::string& outdir, 
        const std::vector<std::string >& references, 
        const std::vector<FramePair >& framepairs, 
        const std::vector<Singlets >& singlets_vec, 
        const std::vector<Multiplets >& multiplets_vec) 
    {   
        SolutionCoder coder;
        for (int32 time = 0; time < framepairs.size(); time ++) {
            LabelAssociations association;
            coder.encode(
                framepairs[time].solution(), framepairs[time].events(),
                singlets_vec[time], singlets_vec[time+1], 
                multiplets_vec[time], multiplets_vec[time+1],
                association);

            // format the output file name
            char file[256];
            sprintf(file, "%08d.txt", time);
/*            std::string file = references[time];
            if (file.find_last_of(".tiff") != std::string::npos)
                file = file.substr(0, file.find_last_of(".tiff"));
            else if (file.find_last_of(".tif") != std::string::npos)
                file = file.substr(0, file.find_last_of(".tif"));
*/
            saveAssociationFile(time, outdir + "/" + std::string(file), association);
        }
    };
};

}

#endif /* __TIFF_READER_WRITER_HXX */
