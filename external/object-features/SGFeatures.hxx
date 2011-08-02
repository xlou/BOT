#ifndef SGFEATURES_HXX
#define SGFEATURES_HXX


#include <vigra/transformimage.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/basicimage.hxx>
#include <vigra/stdimage.hxx>

#include <vigra/functorexpression.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/multi_pointoperators.hxx>

#include <vigra/mathutil.hxx>


#include "ConfigFeatures.hxx"

#include "math.h"
#include "iostream"


#include <vigra/timing.hxx>

namespace SGF {
    using namespace vigra::functor;

    /**
     * Functor to calculate center coordinates in N dimensions.
     * The coordinates are passed as a TinyVector of size N.
     */
    template <class T, int N>
    class FindCenter
    {
    public:

        typedef vigra::TinyVector<T,N> argument_type;
        typedef vigra::TinyVector<double,N> result_type;
        result_type sum_coordinates;
        unsigned int count;

//        typedef typename vigra::MultiArrayView<N,T>::iterator iterator_type;
        typedef typename vigra::MultiArray<N,T>::iterator iterator_type;

        FindCenter()
            : sum_coordinates((double)0), count(0)
        {}

        void reset()
        {
            count = 0;
            sum_coordinates = result_type(0);
        }

        inline void operator()(argument_type const & coord)
        {
            sum_coordinates += coord;
            count++;
        }

        inline void operator()(iterator_type const & it)
        {
            sum_coordinates += it.point();
            count++;
        }

        inline void operator()(FindCenter const & other)
        {
            count += other.count;
            sum_coordinates += other.sum_coordinates;
        }

        inline result_type operator()() const
        {
            if(count != 0)
                return result_type(sum_coordinates / count);
            else
                return result_type(-1);
        }
    };

    /** Masked version of FindCenter
      * Only updates, if the mask is != 0
      */
    template <class T, class MASKTYPE, int N>
    class FindCenterMask: public FindCenter<T,N>
    {
    public:
        typedef MASKTYPE mask_type;
        typedef vigra::TinyVector<T,N> argument_type;
        typedef vigra::TinyVector<T,N> result_type;

//        typedef typename vigra::MultiArrayView<N,T>::iterator iterator_type;
        typedef typename vigra::MultiArray<N,T>::iterator iterator_type;

        inline void operator()(argument_type const & coord, mask_type const & m)
        {
            if((unsigned int)m != 0){
                FindCenter<T,N>::operator()(coord);
            }
        }

        inline void operator()(iterator_type const & it, mask_type const & m)
        {
            if((unsigned int)m != 0){
                FindCenter<T,N>::operator()(it.point());
            }
        }

        inline result_type operator()()
        {
            return FindCenter<T,N>::operator ()();
        }
    };




    /**
     * Functor to calculate the irregularity of a region.
     * For details see the original paper
     */
    template <class T, int N>
    class FindIRGL
    {
    public:

        typedef vigra::TinyVector<T,N> argument_type;
        typedef double result_type;
        result_type IRGL;
        vigra::TinyVector<double,N>  center;
        unsigned int size;

//        typedef typename vigra::MultiArrayView<N,T>::iterator iterator_type;
        typedef typename vigra::MultiArray<N,T>::iterator iterator_type;

        FindIRGL()
            : IRGL(0), center((double)0), size(0)
        {}

        FindIRGL(vigra::TinyVector<double,N>  center, unsigned int size)
            : IRGL(0), center(center), size(size)
        {}

        void reset()
        {
            IRGL = result_type(0);
            center = argument_type(0);
        }


        inline void operator()(argument_type const & coord)
        {
            if(size != 0){
                // FIXME: this is still specialized to 3D. Needs to be extended to ND
                double temp_result = (1+std::pow(4.*M_PI/3., 1./3.)*(coord - center).magnitude())/std::sqrt(size) -1;

                if(temp_result > IRGL){
                    IRGL = temp_result;
                }
            }
        }

        inline void operator()(iterator_type const & it)
        {
            if(size != 0){
                // FIXME: this is still specialized to 3D. Needs to be extended to ND
                double temp_result = (1+std::pow(4.*M_PI/3., 1./3.)*(it.point() - center).magnitude())/std::sqrt(size) -1;

                if(temp_result > IRGL){
                    IRGL = temp_result;
                }
            }
        }

        inline void operator()(FindIRGL const & other)
        {
            center = other.center;
            size = other.size;
            IRGL = other.IRGL;
        }

        inline result_type operator()() const
        {
            return IRGL;
        }
    };



    /** Masked version of vigra::FindMinMax
      * Only updates, if the mask is != 0
      */
    template <class VALUETYPE, class MASKTYPE>
    class FindMinMaxMask: public vigra::FindMinMax<VALUETYPE>
    {
    public:
        typedef MASKTYPE mask_type;
        typedef VALUETYPE argument_type;

        inline void operator()(argument_type const & v, mask_type const & m)
        {
            if((unsigned int)m != 0){
                vigra::FindMinMax<VALUETYPE>::operator()(v);
            }
        }
    };


    /** Inspect a MultiArray and pass coordinates and values to the functor
      *
      * The functor is expected to be able to handle the StridedScanOrderIterator
      * correctly.
    */
    template <int N, class T, class F>
    inline void inspectMultiArrayWithCoordinatesSOI(vigra::MultiArrayView<N,T> &array, F & functor){
//        for(typename vigra::MultiArrayView<N,T>::iterator it = array.begin(); it!=array.end(); ++it )
//            functor( it.point(), *it );
        for (typename vigra::MultiArrayView<N,T >::difference_type_1 ind = 0; ind < array.elementCount(); ind ++) {
            functor(array.scanOrderIndexToCoordinate (ind), array[ind] );
        }
    }

    /** Function to calculate the statistics over all threshold levels
      */
    template <class T>
    feature_array statistics(vigra::MultiArray<1,T> features)
    {
        int len = features.size(0);

        // return values:
        // [0] max value; [1] average value; [2] sample mean; [3] sample stddev
        feature_array ret (array_shape(4),-1.);

        if(len <= 1){
            return ret;
        }

        // calculate statistics by summation
        feature_type max_value  = features[0];
        feature_type sum_feat   = 0.;
        feature_type sum_t_feat = 0.;
        feature_type sum_t_sqr  = 0.;

        for(int t = 0; t < len; t++){
            if(features[t]>max_value)
                max_value = features[t] ;
            sum_feat   += features[t];
            sum_t_feat += (t+1) * features[t];
        }

        ret[0] = max_value;
        ret[1] = sum_feat / (len );

        if(sum_feat == 0){
            return ret;
        }

        ret[2] = sum_t_feat / sum_feat;

        for(int t = 0; t < len; t++){
            sum_t_sqr += std::pow(t+1-ret[2],2) * features[t];
        }

        if(sum_t_sqr < 0){
            return ret;
        }
        ret[3] = std::sqrt( sum_t_sqr / sum_feat );

        return ret;

    }


    template <int N, class T>
    feature_array SGFeatures(vigra::MultiArray<N,T> src_volume, unsigned int levels = 16, T background = 0)
    {

        typedef bool mask_type;
        typedef T value_type;

        // create a mask volume
        vigra::MultiArray<N, mask_type> mask_volume (src_volume.shape());
        vigra::transformMultiArray(vigra::srcMultiArrayRange(src_volume),
                                   vigra::destMultiArray(mask_volume),
                                   vigra::functor::ifThenElse(
                                           Arg1() == Param(background),
                                           Param(0),
                                           Param(1)
                                           )
                                   );




        // get minimum/maximum values
        FindMinMaxMask<value_type,mask_type> minmax;
        vigra::inspectTwoMultiArrays(vigra::srcMultiArrayRange(src_volume),
                                     vigra::srcMultiArray(mask_volume),
                                     minmax
                                     );

        value_type maxint = minmax.max;
        value_type minint = minmax.min;
        value_type intrange = maxint - minint;

        // Storage for the features
        feature_array NCA[2]     = {feature_array(array_shape(levels)),
                                    feature_array(array_shape(levels))};
        feature_array IRGL[2]    = {feature_array(array_shape(levels)),
                                    feature_array(array_shape(levels))};
        feature_array DISP[2]    = {feature_array(array_shape(levels)),
                                    feature_array(array_shape(levels))};
        feature_array INERTIA[2] = {feature_array(array_shape(levels)),
                                    feature_array(array_shape(levels))};
        feature_array TAREA[2]   = {feature_array(array_shape(levels)),
                                    feature_array(array_shape(levels))};
        feature_array CAREA[2]   = {feature_array(array_shape(levels)),
                                    feature_array(array_shape(levels))};

        // apply all thresholds
        for (int t = 0; t < levels; t++)
        {
            label_volume thr_volume (volume_shape(src_volume.shape()));

            double threshold = minint + double(intrange) / double(levels) * (t+1);

            // thresholding
            vigra::transformMultiArray(vigra::srcMultiArrayRange(src_volume),
                                       vigra::destMultiArray(thr_volume),
                                       vigra::functor::ifThenElse(
                                               Arg1() < Param(threshold),
                                               Param(0),
                                               Param(1)
                                               )
                                       );


            unsigned int label_count[2];
            label_volume lbl_volume [2] =
            {label_volume(volume_shape(src_volume.shape())),
             label_volume(volume_shape(src_volume.shape()))};

            // FIXME: this is still specialized to 3D. Needs to be extended to ND
            if(N!=3){
                std::cout << "This is " << N << "D data, I can only handle 3D. So I will crash soon..." << std::endl;
            }
            // first label only the 1-regions (higher than threshold)
            label_count[1] = vigra::labelVolumeWithBackground(
                    vigra::srcMultiArrayRange(thr_volume),
                    vigra::destMultiArray(lbl_volume[1]),
                    vigra::NeighborCode3DSix(),0);

            // relabel all _masked_ threshold labels: 0-->1, 1-->0
            // Idea: Add the mask image to the label image. Areas with values != 1 are mapped to 0
            label_volume sum_volume (src_volume.shape());

            vigra::combineTwoMultiArrays(vigra::srcMultiArrayRange(lbl_volume[1]),
                                       vigra::srcMultiArray(mask_volume),
                                       vigra::destMultiArray(sum_volume),
                                       Arg1()+Arg2() );

            vigra::transformMultiArray(vigra::srcMultiArrayRange(sum_volume),
                                       vigra::destMultiArray(thr_volume),
                                       vigra::functor::ifThenElse(
                                               Arg1() == Param(1),
                                               Param(1),
                                               Param(0) ) );


            // FIXME: this is still specialized to 3D. Needs to be extended to ND
            // label 1-regions (now lower than threshold)
            label_count[0] = vigra::labelVolumeWithBackground(
                    vigra::srcMultiArrayRange(thr_volume),
                    vigra::destMultiArray(lbl_volume[0]),
                    vigra::NeighborCode3DSix(),0);

            // Get the Center Of Gravity for each region
            typedef vigra::ArrayOfRegionStatistics< FindCenter<label_type,N> > f_center;
            f_center COG[2] = {f_center(label_count[0]), f_center(label_count[1])};

            // Get the size of each region
            typedef vigra::ArrayOfRegionStatistics<
                    vigra::FindROISize<label_type> > f_size;
            f_size size[2] = {f_size(label_count[0]),f_size(label_count[1])};

            // Get the irregularity for each region
            typedef vigra::ArrayOfRegionStatistics< FindIRGL<label_type,N> > f_irgl;
            f_irgl IRGL_j[2] = {f_irgl(label_count[0]), f_irgl(label_count[1])};

            // Go through volume and collect the above values
            for (int i = 0; i < 2; ++i)
            {
                inspectMultiArrayWithCoordinatesSOI<N,T,f_center>(lbl_volume[i], COG[i]);

                vigra::inspectTwoMultiArrays(srcMultiArrayRange(lbl_volume[i]),
                                             srcMultiArray(lbl_volume[i]), size[i]);

                // Init irregularity functor with center and size of each clump
                for(int j = 1; j <= label_count[i]; j++){
                    IRGL_j[i][j]( FindIRGL<label_type,N> (COG[i][j](), size[i][j]()) );
                }
                inspectMultiArrayWithCoordinatesSOI<N,T,f_irgl>(lbl_volume[i], IRGL_j[i]);
            }

            // Calculate the Center Of Gravity of the Nucleus
            FindCenterMask<label_type,mask_type,N> COGN;
            inspectMultiArrayWithCoordinatesSOI<N,mask_type,FindCenterMask<label_type,mask_type,N> >(mask_volume,COGN);

            double total_size = COGN.count;
            double sqrt_total_size = std::sqrt(total_size);
            double sqrt_pi = std::sqrt(M_PI);


            // Now calculate the actual feature values:
            //
            // NCA -- Normalised Number of Connected Regions
            // IRGL -- Irregularity
            // DISP -- Average Clump Displacement
            // INERTIA -- Average Clump Interia
            // TAREA -- Total Clump Area
            // CAREA -- Average Clump Area
            for (int i = 0; i < 2; ++i)
            {
                double sum_Dj = 0.;       // normalized clump displacement
                double sum_Dj_NOPj = 0.;  // D_j * NOP_j for average Inertia
                double sum_NOPj = 0.;     // area calculations
                double sum_IRGLj = 0.;    // irregularity

                for (int j = 1; j <= label_count[i]; ++j)
                {
                    double Dj = ( sqrt_pi * (COG[i][j]() - COGN()).magnitude() / sqrt_total_size);
                    sum_Dj += Dj;
                    sum_Dj_NOPj += Dj * size[i][j]();
                    sum_NOPj += size[i][j]();
                    sum_IRGLj += IRGL_j[i][j]()*size[i][j]();
                }

                NCA[i][t] = label_count[i] / total_size;
                IRGL[i][t] = sum_IRGLj/sum_NOPj;
                TAREA[i][t] = sum_NOPj / total_size;

                if(label_count[i] > 0){
                    DISP[i][t] = sum_Dj / label_count[i];
                    INERTIA[i][t] = sum_Dj_NOPj / label_count[i];
                    CAREA[i][t] = sum_NOPj / label_count[i];
                }
            }

        } // END loop over t

        // calculate statistics over all threshold levels
        feature_array STAT_NCA[2]     = {feature_array(array_shape(4)),
                                         feature_array(array_shape(4))};
        feature_array STAT_IRGL[2]    = {feature_array(array_shape(4)),
                                         feature_array(array_shape(4))};
        feature_array STAT_DISP[2]    = {feature_array(array_shape(4)),
                                         feature_array(array_shape(4))};
        feature_array STAT_INERTIA[2] = {feature_array(array_shape(4)),
                                         feature_array(array_shape(4))};
        feature_array STAT_TAREA[2]   = {feature_array(array_shape(4)),
                                         feature_array(array_shape(4))};
        feature_array STAT_CAREA[2]   = {feature_array(array_shape(4)),
                                         feature_array(array_shape(4))};

        for(int i = 0;  i < 2; i++){
            STAT_NCA[i]    = statistics(NCA[i]);
            STAT_IRGL[i]   = statistics(IRGL[i]);
            STAT_DISP[i]   = statistics(DISP[i]);
            STAT_INERTIA[i]= statistics(INERTIA[i]);
            STAT_TAREA[i]  = statistics(TAREA[i]);
            STAT_CAREA[i]  = statistics(CAREA[i]);
        }

        // fill values into output vector
        feature_array ret (array_shape(48));
        for(int i = 0; i < 2; i++){
            int o = i*24;
            for(int j = 0; j < 4; j++){
                ret[o + j]      = STAT_NCA[i][j];
                ret[o + 4 + j]  = STAT_IRGL[i][j];
                ret[o + 8 + j]  = STAT_DISP[i][j];
                ret[o + 12 + j] = STAT_INERTIA[i][j];
                ret[o + 16 + j] = STAT_TAREA[i][j];
                ret[o + 20 + j] = STAT_CAREA[i][j];
            }
        }

        return ret;
    }

}; // END namespace SGF



#endif // SGFEATURES_HXX
