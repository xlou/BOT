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

#ifndef __HDF5_READER_WRITER_HXX
#define __HDF5_READER_WRITER_HXX

#include "vigra/hdf5impex.hxx"
#include "TypeDefinition.hxx"
#include "SolutionCoder.hxx"
#include "TrainingData.hxx"

using namespace vigra::linalg;

namespace bot 
{

/*! This class allows for reading the images/segmentations from and writing 
 *  the tracking results to a hdf5 file.
 *
 */
class HDF5ReaderWriter
{
public:
    /*! A static function that loads the raw images and segmentations
     *  @param filename The path to the hdf5 file
     *  @param images The raw images
     *  @param segmentations The segmentations
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool load(const std::string& filename, 
                   std::vector<Matrix2D >& images,
                   std::vector<Matrix2D >& segmentations) 
    {
        // load raw images
        vigra::HDF5File file(filename, vigra::HDF5File::Open);
        try {
            file.cd("/Raw");
            std::vector<std::string > dirs = file.ls();
            for (int32 time = 0; time < dirs.size(); time ++) {
                file.cd(dirs[time]);
                vigra::ArrayVector<hsize_t > shape = file.getDatasetShape("Data");
                vigra::MultiArray<2, int32 > mat(Shape2D(shape[0], shape[1]), static_cast<int32 >(0));
                file.read<2, int32 >("Data", mat);
                file.cd_up();
                images.push_back(Matrix2D(mat));
            }
        }
        catch (...) {
            std::cerr << "*Warning* Error occured when loading data the raw images" << std::endl;
            return false;
        }
        
        // load segmentations
        try {
            file.cd("/Segmentation");
            std::vector<std::string > dirs = file.ls();
            for (int32 time = 0; time < dirs.size(); time ++) {
                file.cd(dirs[time]);
                vigra::ArrayVector<hsize_t > shape = file.getDatasetShape("Data");
                vigra::MultiArray<2, int32 > mat(Shape2D(shape[0], shape[1]), static_cast<int32 >(0));
                file.read<2, int32 >("Data", mat);
                file.cd_up();
                segmentations.push_back(Matrix2D(mat));
            }
        }
        catch (...) {
            std::cerr << "*Warning* Error occured when loading data the segmentations" << std::endl;
            return false;
        }

        // check compatibility
        if (segmentations.size() != images.size()) {
            std::cerr << "*Warning* The segmentation and the raw image sequence do not share the same length" << std::endl;
            return false;
        }

        return true;
    };

    /*! A static function that loads the training data
     *  @param filename The path to the hdf5 file
     *  @param training The training data
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool load(const std::string& filename,
                   TrainingData& training) 
    {
        // load training data
        vigra::HDF5File file(filename, vigra::HDF5File::Open);
        try {
            file.cd("/Training");
            std::vector<std::string > dirs = file.ls();
            for (int32 ind = 0; ind < dirs.size(); ind ++) {
                int32 time = atoi(dirs[ind].c_str());
                file.cd(dirs[ind]);
                LabelAssociations associations; 

                std::vector<std::string > names = file.ls();
                for (int32 n = 0; n < names.size(); n ++) {
                    LabelAssociation association;
                    association.name = names[n].substr(0, names[n].size()-1);
                    file.cd(names[n]);

                    {
                        vigra::ArrayVector<hsize_t > shape = file.getDatasetShape("Source/Data");
                        vigra::MultiArray<2, int32 > mat(Shape2D(shape[0], shape[1]), static_cast<int32 >(0));
                        file.read<2, int32 >("Source/Data", mat);
                        association.source = Matrix2D(mat);
                    }
                    {
                        vigra::ArrayVector<hsize_t > shape = file.getDatasetShape("Target/Data");
                        vigra::MultiArray<2, int32 > mat(Shape2D(shape[0], shape[1]), static_cast<int32 >(0));
                        file.read<2, int32 >("Target/Data", mat);
                        association.target = Matrix2D(mat);
                    }

                    associations.push_back(association);
                    file.cd_up();
                }
                file.cd_up();

                training.add(time, associations);
            }
        }
        catch (...) {
            std::cerr << "*Warning* Error occured when loading data the training data" << std::endl;
            return false;
        }

        return true;
    };

    /*! A static function that saves the label associations at a given time
     *  @param filename The path to the hdf5 file
     *  @param time The time
     *  @param associations The label associations
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool save(
        const std::string& filename, 
        const int32 time,
        const LabelAssociations& associations)
    {
        char buf[1024];
        // load training data
        vigra::HDF5File file(filename, vigra::HDF5File::Open);
        try {
            file.cd_mk("/Tracking");
            sprintf(buf, "%08d", time);
            file.cd_mk(buf);
            for (int32 ind = 0; ind < associations.size(); ind ++) {
                file.cd_mk(associations[ind].name);
                file.cd_mk("Source");
                vigra::MultiArray<2, int32 > source(associations[ind].source);
                file.write("Data", source);
                file.cd_up();
                file.cd_mk("Target");
                vigra::MultiArray<2, int32 > target(associations[ind].target);
                file.write("Data", target);
                file.cd_up();
                file.cd_up();
            }
            file.cd_up();
        }
        catch (...) {
            std::cerr << "*Warning* Error occured when saving the tracking results" << std::endl;
            return false;
        }

        return true;
    };

    /*! A static function that saves all the tracking solutions
     *  @param filename The path to the hdf5 file
     *  @param framepairs The vector of framepairs
     *  @param singlets_vec The vector of Singlets objects
     *  @param multiplets_vec The vector of Multiplets objects
     *  @return Return true if no error occurs, false otherwise.
     */
    static bool save(
        const std::string& filename, 
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
            save(filename.c_str(), time, association);
        }
    };
};

}

#endif /* __HDF5_READER_WRITER_HXX */
