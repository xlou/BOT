#ifndef OBJECTFEATURES_HXX
#define OBJECTFEATURES_HXX


#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

#include <vigra/multi_array.hxx>
#include <vigra/matrix.hxx>
#include <vigra/eigensystem.hxx>
#include <vigra/labelvolume.hxx>

#include "ConfigFeatures.hxx"
#include "SGFeatures.hxx"


namespace features {

    /** Merge statistics of two objects, e.g. their intensity distribution statistics
      *
      * Input:   weight (w), mean (m), variance (v), skew (s),
      *          and kurtosis (k) of both objects
      * Output:  1D Array, length 5 with new weight (0), mean (1), variance (2),
      *          skew (3), and kurtosis (4)
      * weight is e.g. the number of voxels.
      */
    inline vigra::MultiArray<1,double> mergeStatistics(
            double w1, double m1, double v1, double s1, double k1,
            double w2, double m2, double v2, double s2, double k2
            ){
        vigra::MultiArray<1,double> ret (vigra::MultiArrayShape<1>::type(5));

        ret[0] = w1 + w2;

        double A1 = w1 / ret[0];
        double A2 = w2 / ret[0];

        // merge the mean for every coordinate
        ret[1] = A1 * m1 + A2 * m2;


        // merge the variance for every coordinate
        ret[2] = A1 * (v1 + std::pow( m1 , 2)) + A2 * (v2 + std::pow( m2 , 2))
                 - std::pow( A1 * m1 + A2 * m2 , 2);


        // merge the skew for every coordinate
        ret[3] = A1 * ( std::pow(v1,3./2.)*s1 + 3*m1*v1 + std::pow(m1,3) )
                 +A2* ( std::pow(v2,3./2.)*s2 + 3*m2*v2 + std::pow(m2,3) )
                 - 3*ret[1]*ret[2] - std::pow(ret[1],3);
        ret[3] /= std::pow(ret[2],3./2.);


        // merge the kurtosis for every coordinate
        ret[4] = A1 * ( std::pow(v1,2)*k1 + 4*m1*std::pow(v1,3./2.)*s1 + 6*std::pow(m1,2)*v1 + std::pow(m1,4) )
                 +A2* ( std::pow(v2,2)*k2 + 4*m2*std::pow(v2,3./2.)*s2 + 6*std::pow(m2,2)*v2 + std::pow(m2,4) )
                 -4 * ret[1]*std::pow(ret[2],3./2.)*ret[3] - 6 * std::pow(ret[1],2)*ret[2] - std::pow(ret[1],4) ;
        ret[4] /= std::pow(ret[2],2);

        return ret;

    }




    typedef vigra::MultiArray<1,double> features_type;
    typedef vigra::MultiArray<1,bool> valid_type;


    /** Feature base class
      * All features must be derived from this class.
      */
    template<class T, int N>
    class ObjectFeature{
    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Base class constuctor
          * Only initialize the _length member variable
          */
        ObjectFeature ( int length ) : _length(length) {}


        /** Extract function (not implemented, just as a prototype) to calculate the feature values
          */
        void extract(coordinates_type& coordinates, values_type& values);


        /** getter for the calculated feature values
          */
        features_type get () {
            return _feature_values;
        }


        /** setter for the feature values
          */
        void set (features_type & features) {
            if(features.shape(0) == _length){
                _feature_values = features;
            }
        }


    protected:
        /// number of calculated features
        const int _length;

        /// the calculated feature values
        features_type _feature_values;

        /// indicates if the feature contains valid content
        valid_type _feature_valid;
    };




    /** Object volume feature
      * Count the number of voxels the object consists of.
      * Output stucture:
      * size = 1
      * [0] volume
      */
    template<class T, int N>
    class ObjectVolume: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectVolume(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>(1) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the volume from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the volume from an object (iterator interface)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            int count = 0;
            for(; coordinates_begin != coordinates_end; coordinates_begin++, count++){}

            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length));
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);
            _feature_values(0) = double(count);
            _feature_valid(0) = true;
        }


        /** Merge two objects
          * The volume of the merged object is calculated
          */
        bool mergeWith(ObjectVolume & other) {
            _feature_values(0) = _feature_values(0) + other._feature_values(0);
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;

    };




    /** (Unweighted) Mean position feature
      * Calculate the mean position and its higher central moments. Only the
      * shape of the object is used, not the intensities/values.
      * Output stucture:
      * size = 4*N
      * [0 .. N-1] mean coordinates
      * [N .. 2*N-1] variance
      * [2*N .. 3*N-1] skew
      * [3*N .. 4*N-1] kurtosis
      */
    template<class T, int N>
    class ObjectPosition: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectPosition(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>(4*N) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the position features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }
        /** Extract the position features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            int size = coordinates_end - coordinates_begin;
            _volume = size;

            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            //catch cases with empty or small lists
            if(size == 0){
                return;
            }

            if(size == 1){
                for(int i = 0; i < N; i++){
                    _feature_values[i] = (*coordinates_begin)[i];
                    _feature_valid[i] = true;
                }
                return;
            }

            int n = N;
            features_type sum_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type var_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type skew_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type kurt_coords (vigra::MultiArrayShape<1>::type(n), 0.);

            // calculate mean for all N coordinates
            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int j = 0; j < N; j++){
                    sum_coords[j] += (*it)[j];
                }
            }

            for(int i = 0; i < N; i++){
                _feature_values[i] = sum_coords[i]/size;
                _feature_valid[i] = true;
            }


            // then calculate central moments
            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int j = 0; j < N; j++){
                    double delta = (*it)[j] - _feature_values[j];
                    var_coords[j] += std::pow(delta,2);
                    skew_coords[j] += std::pow(delta,3);
                    kurt_coords[j] += std::pow(delta,4);
                }
            }

            for(int i = 0; i < N; i++){
                // variance
                _feature_values[N+i] = var_coords[i]/size;
                _feature_valid[N+i] = true;

                if(var_coords[i] != 0){
                    // skew
                    _feature_values[2*N+i] = skew_coords[i]/double(size)/std::pow(double(_feature_values[N+i]),3./2.);
                    _feature_valid[2*N+i] = true;

                    // kurtosis
                    _feature_values[3*N+i] = kurt_coords[i]/double(size)/std::pow(double(_feature_values[N+i]),2);
                    _feature_valid[3*N+i] = true;
                }
            }
        }


        /** Merge two objects
          * The position features of the merged object are calculated from the
          * features of the original object.
          * Returns false, if merging was not possible
          */
        bool mergeWith(ObjectPosition & other) {
            features_type merged_features (vigra::MultiArrayShape<1>::type(_length),-1.);
            double merged_volume = _volume + other._volume;

            for(int i = 0; i < N; i++){
                /* Explanation:
                   _feature_values[i] is the mean of coordinate i in this object
                   _feature_values[N+i] is the variance of coordinate i in this object
                   _feature_values[2*N+i] is the skew of coordinate i in this object
                   _feature_values[3*N+i] is the kurtosis of coordinate i in this object

                   for the other object accordingly other._feature_values[...]

                   the mean of the merged object is: merged_features[i]
                   var, skew and kurtosis accordingly
                  */

                features_type merged_statistics = mergeStatistics(_volume, _feature_values[i], _feature_values[N+i], _feature_values[2*N+i], _feature_values[3*N+i],
                                                                  other._volume, other._feature_values[i], other._feature_values[N+i], other._feature_values[2*N+i], other._feature_values[3*N+i] );

                merged_features[i] = merged_statistics[1];
                merged_features[N+i] = merged_statistics[2];
                merged_features[2*N+i] = merged_statistics[3];
                merged_features[3*N+i] = merged_statistics[4];

                if(_feature_valid[i] == false || other._feature_valid[i] == false)
                    _feature_valid[i] = false;

                if(_feature_valid[N+i] == false || other._feature_valid[N+i] == false)
                    _feature_valid[N+i] = false;

                if(_feature_valid[2*N+i] == false || other._feature_valid[2*N+i] == false)
                    _feature_valid[2*N+i] = false;

                if(_feature_valid[3*N+i] == false || other._feature_valid[3*N+i] == false)
                    _feature_valid[3*N+i] = false;
            }

            for(int i = 0; i < _length; i++){
                if (_feature_valid[i] == false)
                    return false;
            }

            _volume = merged_volume;
            _feature_values = merged_features;
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;
        
        /// object volume is needed for merging
        double _volume;

    };




    /** Weighted Mean position feature
      * Calculate the mean position and its higher central moments. Each position
      * is weighted with it's intensity value.
      * Output stucture:
      * size = 4*N
      * [0 .. N-1] mean coordinates
      * [N .. 2*N-1] variance
      * [2*N .. 3*N-1] skew
      * [3*N .. 4*N-1] kurtosis
      */
    template<class T, int N>
    class ObjectWeightedPosition: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectWeightedPosition(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>(4*N) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the weighted position features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            //catch cases with empty or small lists
            if(coordinates.size() == 0 || coordinates.size() != values.size()){
                return;
            }

            if(coordinates.size() == 1){
                for(int i = 0; i < N; i++){
                    _feature_values[i] = coordinates[0][i];
                    _feature_valid[i] = true;
                }
                return;
            }

            int n = N;
            double sum_values = 0.;
            features_type sum_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type var_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type skew_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type kurt_coords (vigra::MultiArrayShape<1>::type(n), 0.);

            // calculate mean for all N coordinates
            for(int i = 0; i < coordinates.size(); i++){
                double value = values[i];
                sum_values += value;
                for(int j = 0; j < N; j++){
                    sum_coords[j] += value*coordinates[i][j];
                }
            }

            _volume = sum_values;

            for(int i = 0; i < N; i++){
                _feature_values[i] = sum_coords[i]/sum_values;
                _feature_valid[i] = true;
            }


            // then calculate central moments
            for(int i = 0; i < coordinates.size(); i++){
                double value = values[i];
                for(int j = 0; j < N; j++){
                    double delta = coordinates[i][j] - _feature_values[j];
                    var_coords[j] += value*std::pow(delta,2);
                    skew_coords[j] += value*std::pow(delta,3);
                    kurt_coords[j] += value*std::pow(delta,4);
                }
            }

            for(int i = 0; i < N; i++){
                // variance
                _feature_values[N+i] = var_coords[i]/sum_values;
                _feature_valid[N+i] = true;

                if(var_coords[i] != 0){
                    // skew
                    _feature_values[2*N+i] = skew_coords[i]/sum_values/std::pow(double(_feature_values[N+i]),3./2.);
                    _feature_valid[2*N+i] = true;

                    // kurtosis
                    _feature_values[3*N+i] = kurt_coords[i]/sum_values/std::pow(double(_feature_values[N+i]),2);
                    _feature_valid[3*N+i] = true;
                }
            }
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the weighted position features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            int size = coordinates_end - coordinates_begin;

            //catch cases with empty or small lists
            if(size == 0){
                return;
            }

            if(size == 1){
                for(int i = 0; i < N; i++){
                    _feature_values[i] = (*coordinates_begin)[i];
                    _feature_valid[i] = true;
                }
                return;
            }

            int n = N;
            double sum_values = 0.;
            features_type sum_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type var_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type skew_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            features_type kurt_coords (vigra::MultiArrayShape<1>::type(n), 0.);

            // calculate mean for all N coordinates
            C it_c = coordinates_begin; V it_v = values_begin;
            for( ; it_c != coordinates_end; it_c++, it_v++){
                double value = (*it_v);
                sum_values += value;
                for(int j = 0; j < N; j++){
                    sum_coords[j] += value*(*it_c)[j];
                }
            }

            _volume = sum_values;

            for(int i = 0; i < N; i++){
                _feature_values[i] = sum_coords[i]/sum_values;
                _feature_valid[i] = true;
            }


            // then calculate central moments
            it_c = coordinates_begin; it_v = values_begin;
            for( ; it_c != coordinates_end; it_c++, it_v++){
                double value = (*it_v);
                for(int j = 0; j < N; j++){
                    double delta = (*it_c)[j] - _feature_values[j];
                    var_coords[j] += value*std::pow(delta,2);
                    skew_coords[j] += value*std::pow(delta,3);
                    kurt_coords[j] += value*std::pow(delta,4);
                }
            }

            for(int i = 0; i < N; i++){
                // variance
                _feature_values[N+i] = var_coords[i]/sum_values;
                _feature_valid[N+i] = true;

                if(var_coords[i] != 0){
                    // skew
                    _feature_values[2*N+i] = skew_coords[i]/sum_values/std::pow(double(_feature_values[N+i]),3./2.);
                    _feature_valid[2*N+i] = true;

                    // kurtosis
                    _feature_values[3*N+i] = kurt_coords[i]/sum_values/std::pow(double(_feature_values[N+i]),2);
                    _feature_valid[3*N+i] = true;
                }
            }
        }

        /** Merge two objects
          * The position features of the merged object are calculated from the
          * features of the original object.
          * Returns false, if merging was not possible
          */
        bool mergeWith(ObjectWeightedPosition & other) {
            features_type merged_features (vigra::MultiArrayShape<1>::type(_length),-1.);
            double merged_volume = _volume + other._volume;

            double A1 = double(_volume) / double(merged_volume);
            double A2 = double(other._volume) / double(merged_volume);

            for(int i = 0; i < N; i++){
                /* Explanation:
                   _feature_values[i] is the mean of coordinate i in this object
                   _feature_values[N+i] is the variance of coordinate i in this object
                   _feature_values[2*N+i] is the skew of coordinate i in this object
                   _feature_values[3*N+i] is the kurtosis of coordinate i in this object

                   for the other object accordingly other._feature_values[...]

                   the mean of the merged object is: merged_features[i]
                   var, skew and kurtosis accordingly
                  */

                features_type merged_statistics = mergeStatistics(_volume, _feature_values[i], _feature_values[N+i], _feature_values[2*N+i], _feature_values[3*N+i],
                                                                  other._volume, other._feature_values[i], other._feature_values[N+i], other._feature_values[2*N+i], other._feature_values[3*N+i] );

                merged_features[i] = merged_statistics[1];
                merged_features[N+i] = merged_statistics[2];
                merged_features[2*N+i] = merged_statistics[3];
                merged_features[3*N+i] = merged_statistics[4];

                if(_feature_valid[i] == false || other._feature_valid[i] == false)
                    _feature_valid[i] = false;

                if(_feature_valid[N+i] == false || other._feature_valid[N+i] == false)
                    _feature_valid[N+i] = false;

                if(_feature_valid[2*N+i] == false || other._feature_valid[2*N+i] == false)
                    _feature_valid[2*N+i] = false;

                if(_feature_valid[3*N+i] == false || other._feature_valid[3*N+i] == false)
                    _feature_valid[3*N+i] = false;
            }

            for(int i = 0; i < _length; i++){
                if (_feature_valid[i] == false)
                    return false;
            }

            _volume = merged_volume;
            _feature_values = merged_features;
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;

        /// object volume (=sum over all value weights) is needed for merging
        double _volume;

    };




    /** Principal components feature
      * Calculate the major axes (principal components) of object shape
      * Output stucture:
      * size = N*(N+1)
      * [0 .. N-1] Eigenvalues of covariance matrix
      * [k*(N+1) .. k*(N+2)-1] Eigenvector for eigenvalue [k]
      */
    template<class T, int N>
    class ObjectPrincipalComponents: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectPrincipalComponents(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>( N*(N+1) ) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the PC features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            //catch cases with empty or small lists
            if(coordinates.size() <= 1 || coordinates.size() != values.size()){
                return;
            }

            int size = coordinates.size();

            // calculate mean for all N coordinates
            int n = N;
            features_type mean_coords (vigra::MultiArrayShape<1>::type(n), 0.);

            for(int i = 0; i < coordinates.size(); i++){
                for(int j = 0; j < N; j++){
                    mean_coords[j] += coordinates[i][j]/double(size);
                }
            }

            // set up vigra matrix with coordinates shifted to zero mean
            vigra::linalg::Matrix<double> coord_matrix (size,N,0.);
            for(int i = 0; i < size; i++){
                for(int j = 0; j < N; j++){
                    coord_matrix(i,j) = double(coordinates[i][j] - mean_coords[j]);
                }
            }

            // claculate covariance matrix
            vigra::linalg::Matrix<double> cov_matrix (N,N,0.);
            vigra::linalg::covarianceMatrixOfColumns(coord_matrix, cov_matrix);

            // calculate eigenvalues and eigenvectors
            vigra::linalg::Matrix<double> ev_matrix (N,N,0.);	//matrix for eigenvectors
            vigra::linalg::Matrix<double> ew_vector (N,1,0.);	//matrix for eigenvalues

            vigra::linalg::symmetricEigensystem(cov_matrix,ew_vector,ev_matrix);


            // first fill the eigenvalues to the output vector
            for(int i = 0; i < N; i++){
                _feature_values[i] = ew_vector(i,0);
                _feature_valid[i] = true;

                for(int j = 0; j < N; j++){
                    _feature_values[(i+1)*N + j] = ev_matrix(j,i);
                    _feature_valid[(i+1)*N + j] = true;
                }
            }
            extract(coordinates.begin(), coordinates.end(), values.end());
        }

        /** Extract the PC features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            int size = coordinates_end - coordinates_begin;
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            //catch cases with empty or small lists
            if(size <= 1){
                return;
            }

            // calculate mean for all N coordinates
            int n = N;
            features_type mean_coords (vigra::MultiArrayShape<1>::type(n), 0.);

            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int j = 0; j < N; j++){
                    mean_coords[j] += (*it)[j]/double(size);
                }
            }

            // set up vigra matrix with coordinates shifted to zero mean
            vigra::linalg::Matrix<double> coord_matrix (size,N,0.);
            int i = 0; // not beautiful, redo if you want
            for(C it = coordinates_begin; it != coordinates_end; it++, i++){
                for(int j = 0; j < N; j++){
                    coord_matrix(i,j) = double((*it)[j] - mean_coords[j]);
                }
            }

            // claculate covariance matrix
            vigra::linalg::Matrix<double> cov_matrix (N,N,0.);
            vigra::linalg::covarianceMatrixOfColumns(coord_matrix, cov_matrix);

            // calculate eigenvalues and eigenvectors
            vigra::linalg::Matrix<double> ev_matrix (N,N,0.);	//matrix for eigenvectors
            vigra::linalg::Matrix<double> ew_vector (N,1,0.);	//matrix for eigenvalues

            vigra::linalg::symmetricEigensystem(cov_matrix,ew_vector,ev_matrix);


            // first fill the eigenvalues to the output vector
            for(int i = 0; i < N; i++){
                _feature_values[i] = ew_vector(i,0);
                _feature_valid[i] = true;

                for(int j = 0; j < N; j++){
                    _feature_values[(i+1)*N + j] = ev_matrix(j,i);
                    _feature_valid[(i+1)*N + j] = true;
                }
            }
        }

        /** Merge two objects
          * The Principal Components feature needs to be recalculated on the
          * merged object. Thus mergeWith returns always false, the calculated
          * features are not changed.
          */
        bool mergeWith(ObjectPrincipalComponents & other) {
            std::cout << "Warning: Merging not possible for Principal Components feature. Need to be recalculated." << std::endl;
            return false;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;

    };



    /** Bounding Box feature
      * Find the smallest possible box that contains the whole object
      * Output stucture:
      * size = 2*N + 1
      * [0 .. N-1] minimum coordinates (included by the object)
      * [N .. 2*N-1] maximum coordinates (excluded by the object)
      * [2*N] Fill factor: <object volume> / <bounding box volume>
      */
    template<class T, int N>
    class ObjectBoundingBox: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectBoundingBox(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>( 2*N+1 ) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the bounding box features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.end());
        }

        /** Extract the bounding box features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);
            int size = coordinates_end - coordinates_begin;

            //catch cases with empty or small lists
            if(size <= 1){
                return;
            }

            _volume = size;

            // init with some real data value
            for(int i = 0; i < N; i++){
                _feature_values[i] = (*coordinates_begin)[i];	// set min/max to some value
                _feature_values[N+i] = (*coordinates_begin)[i];	// set min/max to some value
            }

            // go through list and find extremal cordinates
            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int i = 0; i < N; i++){
                    if((*it)[i] < _feature_values[i]){
                        _feature_values[i] = (*it)[i]; // update minimum
                    }
                    if((*it)[i] + 1 > _feature_values[N+i]){
                        _feature_values[N+i] = (*it)[i] + 1; // update maximum
                    }
                }
            }

            // calculate fill factor
            _feature_values[2*N] = _volume;
            for(int i = 0; i < N; i++){
                double blength = _feature_values[N+i] - _feature_values[i] + 1;
                if (blength != 0) {
                    _feature_values[2*N] /= blength;
                }else{
                    _feature_valid[2*N] = false;
                    return;
                }
            }
            _feature_valid[2*N] = true;
        }


        /** Merge two objects
          */
        bool mergeWith(ObjectBoundingBox & other) {
            features_type new_values (vigra::MultiArrayShape<1>::type(_length),-1.);
            for(int i = 0; i < N; i++){
                if(_feature_valid[i] == false || other._feature_valid[i] == false || _feature_valid[N+i] == false || other._feature_valid[N+i] == false){
                    return false;
                }
                new_values[i] = std::min(_feature_values[i], other._feature_values[i]);
                new_values[N+i] = std::max(_feature_values[N+i], other._feature_values[N+i]);
            }

            // calculate fill factor
            new_values[2*N] = _volume + other._volume;
            for(int i = 0; i < N; i++){
                double blength = new_values[N+i] - new_values[i] + 1;
                if (blength != 0) {
                    new_values[2*N] /= blength;
                }else{
                    _feature_valid[2*N] = false;
                }
            }

            _feature_values = new_values;
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;

        double _volume;
    };
    
    
    
    
    /** Intensity feature
      * Calculate the mean intensity and its central moments of all object
      * Output stucture:
      * size = 5
      * [0] Mean of intensity distribution
      * [1] Variance of intensity distribution
      * [2] Skew of Intensity distribution
      * [3] Kurtosis of Intensity distribution
      * [4] Kurtosis of Intensity distribution
      */
    template<class T, int N>
    class ObjectIntensity: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectIntensity(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>(5) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the intensity features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the intensity features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            int size = coordinates_end - coordinates_begin;
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            //catch cases with empty or small lists
            if(size == 0){
                return;
            }

            if(size == 1){
                _feature_values[0] = (*values_begin);
                _feature_valid[0] = true;
                return;
            }

            _volume = size;

            //calculate mean first
            double sum_values = 0.;

            C it_c = coordinates_begin;
            V it_v = values_begin;
            for(; it_c != coordinates_end; it_c++ ,it_v++){
                sum_values += double((*it_v));
            }

            _feature_values[0] = sum_values / double(size);
            _feature_valid[0] = true;
            
            _feature_values[4] = sum_values;
            _feature_valid[4] = true;


            // Calculate variance, skew and kurtosis
            double var = 0;
            double skew = 0;
            double kurt = 0;

            it_c = coordinates_begin;
            it_v = values_begin;
            for(; it_c != coordinates_end; it_c++ ,it_v++){
                double delta = double(*it_v) - _feature_values[0];
                var += std::pow(delta,2);
                skew += std::pow(delta,3);
                kurt += std::pow(delta,4);
            }

            _feature_values[1] = var / double(size);
            _feature_valid[1] = true;

            if(_feature_valid[1] != 0){
                _feature_values[2] = skew/double(size)/std::pow(_feature_values[1],3./2.);
                _feature_values[3] = kurt/double(size)/std::pow(_feature_values[1],2);
                _feature_valid[2] = true;
                _feature_valid[3] = true;
            }

        }


        /** Merge two objects
          * The intensity features of the merged object are calculated from the
          * features of the original object.
          * Returns false, if merging was not possible
          */
        bool mergeWith(ObjectIntensity & other) {
            features_type merged_features (vigra::MultiArrayShape<1>::type(_length),-1.);
            double merged_volume = _volume + other._volume;

            features_type merged_statistics = mergeStatistics(_volume, _feature_values[0], _feature_values[1], _feature_values[2], _feature_values[3],
                                                              other._volume, other._feature_values[0], other._feature_values[1], other._feature_values[2], other._feature_values[3] );

            merged_features[0] = merged_statistics[1];
            merged_features[1] = merged_statistics[2];
            merged_features[2] = merged_statistics[3];
            merged_features[3] = merged_statistics[4];

            for(int i = 0; i < _length; i++){
                if(_feature_valid[i] == false || other._feature_valid[i] == false)
                    _feature_valid[i] = false;
            }

            for(int i = 0; i < _length; i++){
                if (_feature_valid[i] == false)
                    return false;
            }

            _volume = merged_volume;
            _feature_values = merged_features;
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;

        /// object volume (=sum over all value weights) is needed for merging
        double _volume;

    };


    /** Minimum/Maximum Intensity feature
      * Find the minimum and the maximum intensity of a object and find the
      * quantiles of the intensity distribution.
      * Output stucture:
      * size = 9
      * [0] Minimum intensity
      * [1] Maximum intensity
      * [2] 5% quantile
      * [3] 10% quantile
      * [4] 25% quantile
      * [5] 50% quantile
      * [6] 75% quantile
      * [7] 90% quantile
      * [8] 95% quantile
      */
    template<class T, int N>
    class ObjectMinMaxIntensity: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectMinMaxIntensity(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>( 9 ) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the intensity min/max/quantiles features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the intensity min/max/quantiles features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            int size = coordinates_end - coordinates_begin;

            //catch cases with empty or small lists
            if(size <= 1){
                return;
            }

            std::sort(values_begin,values_begin+size);

            // minimum/maximum value
            _feature_values[0] = *values_begin;
            _feature_valid[0] = true;
            _feature_values[1] = *(values_begin+size-1) ;
            _feature_valid[1] = true;

            // find the 5%, 10%, 25%, 50%, 75%, 90%, 95% quanitles
            int i = 0;
            V it = values_begin;

            // 5% quantile
            while(i < 0.05*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[2] = *it;
            _feature_valid[2] = true;

            // 10% quantile
            while(i < 0.1*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[3] = *it;
            _feature_valid[3] = true;

            // 25% quantile
            while(i < 0.25*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[4] = *it;
            _feature_valid[4] = true;

            // 50% quantile
            while(i < 0.5*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[5] = *it;
            _feature_valid[5] = true;

            // 75% quantile
            while(i < 0.75*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[6] = *it;
            _feature_valid[6] = true;

            // 90% quantile
            while(i < 0.9*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[7] = *it;
            _feature_valid[7] = true;

            // 95% quantile
            while(i < 0.95*size && i < size-1){
                i++;
                it++;
            }

            _feature_values[8] = *it;
            _feature_valid[8] = true;

        }


        /** Merge two objects
          * Merging is not (yet) implemented for the MinMaxIntensity feature.
          */
        bool mergeWith(ObjectMinMaxIntensity & other) {
            std::cout << "Warning: Merging for MinMaxIntensity Feature is not available." << std::endl;
            return false;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;
    };




    /** Minimum/Maximum Intensity feature
      * Find the minimum and the maximum intensity of a object and find the
      * quantiles of the intensity distribution.
      * Output stucture:
      * size = N+1
      * [0] Maximum intensity
      * [1 .. N] coordinates of max. intensity
      */
    template<class T, int N>
    class ObjectMaxIntensity: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectMaxIntensity(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>( N+1 ) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the MaxIntensity features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the intensity max intensity position features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);
            int size = coordinates_end - coordinates_begin;

            //catch cases with empty lists
            if(size < 1){
                return;
            }

            // calculate mean for all N coordinates
            int n = N;
            features_type mean_coords (vigra::MultiArrayShape<1>::type(n), 0.);

            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int j = 0; j < N; j++){
                    mean_coords[j] += (*it)[j]/double(size);
                }
            }

            // some start values
            features_type max_coords (vigra::MultiArrayShape<1>::type(n), 0.);
            for(int i = 0; i < N; i++){
                max_coords[i] = (*coordinates_begin)[i];
            }
            double max = *values_begin;


            for(;coordinates_begin != coordinates_end; coordinates_begin++, values_begin++)
            {
                if(*values_begin > max){

                    max = *values_begin;
                    for(int j = 0; j < N; j++){
                        max_coords[j] = (*coordinates_begin)[j];
                    }

                }else if(*values_begin == max){

                    double d_old = 0.;
                    double d_new = 0.;
                    for(int j = 0; j < N; j++){
                        d_old += std::pow(max_coords[j]-mean_coords[j],2);
                        d_new += std::pow((*coordinates_begin)[j]-mean_coords[j],2);
                    }

                    d_old = std::sqrt(d_old);
                    d_new = std::sqrt(d_new);

                    if(d_new < d_old){
                        max = *values_begin;
                        for(int j = 0; j < N; j++){
                            max_coords[j] = (*coordinates_begin)[j];
                        }
                    }
                }
            }

            _feature_values[0] = max;
            for(int j = 0; j < N; j++){
                _feature_values[j+1] = max_coords[j];
            }

            _feature_valid[0] = true;
            _feature_valid[1] = true;
            _feature_valid[2] = true;
            _feature_valid[3] = true;
        }


        /** Merge two objects
          */
        bool mergeWith(ObjectMaxIntensity & other) {
            features_type merged_features (vigra::MultiArrayShape<1>::type(_length),-1.);

            if(_feature_values[0] > other._feature_values[0]){
                for(int i = 0; i < _length; i++){
                    merged_features[i] = _feature_values[i];
                }
            }else{
                for(int i = 0; i < _length; i++){
                    merged_features[i] = other._feature_values[i];
                }
            }

            for(int i = 0; i < _length; i++){
                if(_feature_valid[i] == false || other._feature_valid[i] == false)
                    _feature_valid[i] = false;
            }

            for(int i = 0; i < _length; i++){
                if (_feature_valid[i] == false)
                    return false;
            }

            _feature_values = merged_features;
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;
    };



    /** Pairwise features
      * Calculate average values of differences of neighboring intensity values.
      * Attention: Features are not fixed yet.
      * Output stucture:
      * size = 4
      * [0] Average sum over absolute distances
      * [1] Average sum over squared distances
      * [2] Average symmetric first derivative
      * [3] Average second derivative
      */
    // Pairwise features are deprecated. First calculate the absolute/squared/
    // symmetric first derivative and second derivative on the whole volume and
    // then use the ObjectStatistics Feature to calculate statistics over the
    // object.
    template<class T, int N>
    class ObjectPairwise: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectPairwise(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>( 4 ) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the Pairwise features from an object, works only in 3D
          */
        void extract(coordinates_type& coordinates, values_type& values) {
            extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the Pairwise features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),0.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);
            int size = coordinates_end - coordinates_begin;

            //catch cases with empty lists
            if(size < 1){
                return;
            }



            //fill a zero-padded volume
            typename vigra::MultiArrayShape<2*N>::type minmax;
            // init with some real data value
            for(int i = 0; i < 2*N; i++){
                minmax[i] = (*coordinates_begin)[i];	// set min/max to some value
            }

            // go through list and find extremal coordinates
            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int i = 0; i < N; i++){
                    if((*it)[i] < minmax[i]){
                        minmax[i] = (*it)[i]; // update minimum
                    }
                    if((*it)[i] > minmax[N+i]){
                        minmax[N+i] = (*it)[i]; // update maximum
                    }
                }
            }

            typename vigra::MultiArrayShape<N>::type shape;
            for(int i = 0; i < N; i++){
                shape[i] = minmax[N+i] - minmax[i] +1;
            }

            vigra::MultiArray<N,T> volume (shape,T(0));

            C it_c = coordinates_begin;
            V it_v = values_begin;
            for(; it_c != coordinates_end; it_c++, it_v++)
            {
                typename vigra::MultiArrayShape<N>::type pos;
                for(int j = 0; j < N; j++){
                    pos[j] = (*it_c)[j] - minmax[j];
                }

                volume[pos] = *it_v;
            }


            // FIXME: Code is only for 3D from here on. Needs to be implemented for ND or at least 2D
            vigra_precondition(N==3, "Error: Pairwise Features are only implemented in 3D.");
            //go over volume and sum up pairwise costs
            for(int i = 0; i < volume.size(0); i++){
                for(int j = 0; j < volume.size(1); j++){
                    for(int k = 0; k < volume.size(2); k++){

                        if(i+1 < volume.size(0)){
                            // relative sum over absolute distances
                            _feature_values[0] += std::abs(volume(i+1,j,k) - volume(i,j,k))/double(size);
                            // relative sum over squared distances
                            _feature_values[1] += std::pow(volume(i+1,j,k) - volume(i,j,k),2)/double(size);
                            if(i>0){
                                // sum over symmetric derivatives
                                _feature_values[2] += std::abs(volume(i+1,j,k) - volume(i-1,j,k))/double(size);
                                // second derivative
                                _feature_values[3] += std::abs(2*volume(i,j,k) - volume(i-1,j,k) - volume(i+1,j,k))/double(size);
                            }
                        }
                        if(j+1 < volume.size(1)){
                            _feature_values[0] += std::abs(volume(i,j+1,k) - volume(i,j,k))/double(size);
                            _feature_values[1] += std::pow(volume(i,j+1,k) - volume(i,j,k),2)/double(size);
                            if(j>0){
                                _feature_values[2] += std::abs(volume(i,j+1,k) - volume(i,j-1,k))/double(size);
                                _feature_values[3] += std::abs(2*volume(i,j,k) - volume(i,j-1,k) - volume(i,j+1,k))/double(size);
                            }
                        }
                        if(k+1 < volume.size(2)){
                            _feature_values[0] += std::abs(volume(i,j,k+1) - volume(i,j,k))/double(size);
                            _feature_values[1] += std::pow(volume(i,j,k+1) - volume(i,j,k),2)/double(size);
                            if(k>0){
                                _feature_values[2] += std::abs(volume(i,j,k+1) - volume(i,j,k-1))/double(size);
                                _feature_values[3] += std::abs(2*volume(i,j,k) - volume(i,j,k-1) - volume(i,j,k+1))/double(size);
                            }
                        }
                    }
                }
                _feature_valid[0] = true;
                _feature_valid[1] = true;
                _feature_valid[2] = true;
                _feature_valid[3] = true;
            }
        }

        /** Merge two objects
          * Not (yet) implemented for Pairwise features
          */
        bool mergeWith(ObjectPairwise & other) {
            return false;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;
    };




    /** SGF Features
      * Statistical Geometric Features as described in
      * Jackway, Walker: Statistical Geometric Features - Extensions For
      * Cytological Texture Analysis (ICPR, 1996).
      *
      * Output stucture:
      * size = 48
      * The detailed information about the meaning of each field can be found in
      * the original paper.
      */
    template<class T, int N>
    class ObjectSGF: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<typename vigra::MultiArrayShape<N>::type> coordinates_type;
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectSGF(coordinates_type& coordinates, values_type& values): ObjectFeature<T,N>( 48 ) {
            extract(coordinates, values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the SGF features from an object
          */
        void extract(coordinates_type& coordinates, values_type& values) {
             extract(coordinates.begin(), coordinates.end(), values.begin());
        }

        /** Extract the SGF features from an object (iterator version)
          */
        template <class C, class V> void extract(C coordinates_begin, C coordinates_end, V values_begin) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);
            int size = coordinates_end - coordinates_begin;

            //catch cases with empty lists
            if(size < 1){
                return;
            }



            //fill a zero-padded volume
            typename vigra::MultiArrayShape<2*N>::type minmax;
            // init with some real data value
            for(int i = 0; i < N; i++){
                minmax[i] = (*coordinates_begin)[i];	// set min/max to some value
                minmax[N+i] = (*coordinates_begin)[i];
            }

            // go through list and find extremal coordinates
            for(C it = coordinates_begin; it != coordinates_end; it++){
                for(int i = 0; i < N; i++){
                    if((*it)[i] < minmax[i]){
                        minmax[i] = (*it)[i]; // update minimum
                    }
                    if((*it)[i] > minmax[N+i]){
                        minmax[N+i] = (*it)[i]; // update maximum
                    }
                }
            }

            typename vigra::MultiArrayShape<N>::type shape;
            for(int i = 0; i < N; i++){
                shape[i] = minmax[N+i] - minmax[i] +1;
            }

            vigra::MultiArray<N,T> volume (shape,T(0));

            C it_c = coordinates_begin;
            V it_v = values_begin;
            for(; it_c != coordinates_end; it_c++, it_v++)
            {
                typename vigra::MultiArrayShape<N>::type pos;
                for(int j = 0; j < N; j++){
                    pos[j] = (*it_c)[j] - minmax[j];
                }

                volume[pos] = *it_v;
            }


            // FIXME: Code is only for 3D from here on. Needs to be implemented for ND or at least 2D
            vigra_precondition(N==3, "Error: Statistical Geometric Features only work in 3D. A ND connected components implementation is needed.");
            if(N==3){
                // Calculate the Features
                _feature_values = SGF::SGFeatures<N,T>(volume, 24);
            }
        }


        /** Merge two objects
          * Not (yet) implemented for statistical geometric features
          */
        bool mergeWith(ObjectSGF & other) {
            return false;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;
    };



    /** General N-D Statistics feature
      *
      * This feature calculates mean, variance, skew, kurtosis, min, and max of
      * the value distribution. This feature allows merging. No information about
      * the spatial relation is incorporated, thus the interface is simplified.
      */
    template<class T, int N>
    class ObjectStatistics: public ObjectFeature<T,N> {

    public:
        typedef vigra::ArrayVector<T> values_type;


        /** Constructor to calculate feature immediately
          */
        ObjectStatistics(values_type& values): ObjectFeature<T,N>(6) {
            extract(values);
        }


        using ObjectFeature<T,N>::set;
        using ObjectFeature<T,N>::get;


        /** Extract the intensity features from an object
          */
        void extract(values_type& values) {
            _feature_values = features_type(vigra::MultiArrayShape<1>::type(_length),-1.);
            _feature_valid = valid_type(vigra::MultiArrayShape<1>::type(_length),false);

            //catch cases with empty or small lists
            if(values.size() == 0 ){
                return;
            }

            _feature_values[4] = values[0]; // min
            _feature_values[5] = values[0]; // max
            _feature_valid[4] = true;
            _feature_valid[5] = true;

            if(values.size() == 1){
                _feature_values[0] = values[0]; // mean
                _feature_valid[0] = true;
                return;
            }

            int size = values.size();
            _volume = size;

            //calculate mean first
            double sum_values = 0.;

            for(int i = 0; i < size; i++){
                double v = values[i];
                sum_values += v;

                // get minimum/maximum on the way
                if(v < _feature_values[4])
                    _feature_values[4] = v;
                if(v > _feature_values[5])
                    _feature_values[5] = v;
            }

            _feature_values[0] = sum_values / double(size);
            _feature_valid[0] = true;

            // Calculate variance, skew and kurtosis
            double var = 0;
            double skew = 0;
            double kurt = 0;

            for(int i = 0; i < size; i++){
                double delta = double(values[i]) - _feature_values[0];
                var += std::pow(delta,2);
                skew += std::pow(delta,3);
                kurt += std::pow(delta,4);
            }

            _feature_values[1] = var / double(size);
            _feature_valid[1] = true;

            if(_feature_valid[1] != 0){
                _feature_values[2] = skew/double(size)/std::pow(_feature_values[1],3./2.);
                _feature_values[3] = kurt/double(size)/std::pow(_feature_values[1],2);
                _feature_valid[2] = true;
                _feature_valid[3] = true;
            }
        }


        /** Merge two objects
          * The statistics do not need to be recalculated.
          */
        bool mergeWith(ObjectStatistics & other) {
            features_type merged_features (vigra::MultiArrayShape<1>::type(_length),-1.);
            double merged_volume = _volume + other._volume;

            features_type merged_statistics = mergeStatistics(_volume, _feature_values[0], _feature_values[1], _feature_values[2], _feature_values[3],
                                                              other._volume, other._feature_values[0], other._feature_values[1], other._feature_values[2], other._feature_values[3] );

            merged_features[0] = merged_statistics[1];
            merged_features[1] = merged_statistics[2];
            merged_features[2] = merged_statistics[3];
            merged_features[3] = merged_statistics[4];

            merged_features[4] = std::min(_feature_values[4], other._feature_values[4]);
            merged_features[5] = std::max(_feature_values[5], other._feature_values[5]);

            for(int i = 0; i < _length; i++){
                if(_feature_valid[i] == false || other._feature_valid[i] == false)
                    _feature_valid[i] = false;
            }

            for(int i = 0; i < _length; i++){
                if (_feature_valid[i] == false)
                    return false;
            }

            _volume = merged_volume;
            _feature_values = merged_features;
            return true;
        }


    private:
        using ObjectFeature<T,N>::_length;
        using ObjectFeature<T,N>::_feature_values;
        using ObjectFeature<T,N>::_feature_valid;

        /// object volume (=sum over all value weights) is needed for merging
        double _volume;

    };






    /** Object volume feature
      * Count the number of voxels the object consists of.
      * Output stucture:
      * size = 1
      * [0] volume
      */
    feature_array extractVolume(three_set& coordinates, value_set& intensities);


    /** (Unweighted) Mean position feature
      * Calculate the mean position and its higher central moments
      * Output stucture:
      * size = 12
      * [0..2] mean x,y,z coordinates
      * [3..5] variance of x,y,z coordinates
      * [6..8] skew of x,y,z coordinates
      * [9..11] kurtosis of x,y,z coordinates
      */
    feature_array extractPosition(three_set& coordinates, value_set& intensities);


    /** Weighted Mean position feature
      * Calculate the intensity weighted mean position and its higher central moments
      * Output stucture:
      * size = 12
      * [0..2] weighted mean x,y,z coordinates
      * [3..5] variance of x,y,z coordinates
      * [6..8] skew of x,y,z coordinates
      * [9..11] kurtosis of x,y,z coordinates
      */
    feature_array extractWeightedPosition(three_set& coordinates, value_set& intensities);


    /** Principal components feature
      * Calculate the principal components of the voxel distribution
      * Output stucture:
      * size = 12
      * [0..2] Eigenvalues of covariance matrix
      * [3..5] Eigenvector of eigenvalue [0]
      * [6..8] Eigenvector of eigenvalue [1]
      * [9..11] Eigenvector of eigenvalue [2]
      */
    feature_array extractPrincipalComponents(three_set& coordinates, value_set& intensities);


    /** Bounding Box feature
      * Find the smallest possible box that contains the whole object
      * Output stucture:
      * size = 7
      * [0] Lower x position
      * [1] Lower y position
      * [2] Lower z position
      * [3] Upper x position
      * [4] Upper y position
      * [5] Upper z position
      * [6] Fill factor: <object volume> / <size of bounding box>
      */
    feature_array extractBoundingBox(three_set& coordinates, value_set& intensities);


    /** Intensity feature
      * Calculate the mean intensity and its central moments of all object
      * Output stucture:
      * size = 4
      * [0] Mean of intensity distribution
      * [1] Variance of intensity distribution
      * [2] Skew of Intensity distribution
      * [3] Kurtosis of Intensity distribution
      */
    feature_array extractIntensity(three_set& coordinates, value_set& intensities);


    /** Minimum/Maximum Intensity feature
      * Find the minimum and the maximum intensity of a object and find the
      * quantiles of the intensity distribution.
      * Output stucture:
      * size = 9
      * [0] Minimum intensity
      * [1] Maximum intensity
      * [2] 5% quantile
      * [3] 10% quantile
      * [4] 20% quantile
      * [5] 50% quantile
      * [6] 80% quantile
      * [7] 90% quantile
      * [8] 95% quantile
      */
    feature_array extractMinMaxIntensity(three_set& coordinates, value_set& intensities);

    /** Minimum/Maximum Intensity feature
      * Find the minimum and the maximum intensity of a object and find the
      * quantiles of the intensity distribution.
      * Output stucture:
      * size = 4
      * [0] Maximum intensity
      * [1] x1 coord
      * [2] x2 coord
      * [3] x3 coord
      */
    feature_array extractMaxIntensity(three_set& coordinates, value_set& intensities);



    /** Pairwise features
      * Calculate average values of differences of neighboring intensity values.
      * Attention: Features are not fixed yet.
      * Output stucture:
      * size = 4
      * [0] Average sum over absolute distances
      * [1] Average sum over squared distances
      * [2] Average symmetric first derivative
      * [3] Average second derivative
      */
    feature_array extractPairwise(three_set &coordinates, value_set &intensities);


    /** Histogram of Oriented Gradients
      * Don't use them. They don't improve anything, just cost computation time.
      * RF variable importance was about 0.1 of informative features
      * Output stucture:
      * size = 20
      */
    feature_array extractHOG(three_set &coordinates, value_set &intensities);


    /** Experimental SGF Features
      *
      * Output stucture:
      * size = 48
      */
    feature_array extractSGF(three_set &coordinates, value_set &intensities);


} /* namespace features */
#endif // OBJECTFEATURES_HXX
