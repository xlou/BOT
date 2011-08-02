
#ifndef CONFIGFEATURES_HXX
#define CONFIGFEATURES_HXX

#include <vector>

#include <vigra/multi_array.hxx>
#include <vigra/matrix.hxx>
#include <vigra/array_vector.hxx>

// Datatypes used for raw data

typedef unsigned short raw_type;
typedef vigra::MultiArray<3,raw_type> raw_volume;


// Datatypes used for segmentation

typedef unsigned short seg_type;
typedef vigra::MultiArray<3,seg_type> seg_volume;

// Datatypes used for labeled volume

typedef unsigned short label_type;
typedef vigra::MultiArray<3,label_type> label_volume;


// Shape of 3D MultiArray

typedef vigra::MultiArrayShape<3>::type volume_shape;


// Datatype for coordinates

typedef short coordinate_type;

// 2D shape
typedef vigra::MultiArrayShape<2>::type matrix_shape;

// 1D MultiArrays / Lists

typedef vigra::MultiArray<1,label_type> label_array;  // arrays containing labels

typedef vigra::MultiArray<1,coordinate_type> dim_size_type; // arrays containing coordinate values

typedef float feature_type; // all feature values have this type
typedef vigra::MultiArray<1,feature_type> feature_array; // list of features

typedef vigra::MultiArrayShape<1>::type array_shape;

typedef vigra::MultiArray<2, feature_type> feature_matrix;


// Datatype for coordinate sets (vigra-based -> easy to read/write)

typedef vigra::MultiArray<2,coordinate_type> set_type;
typedef vigra::MultiArrayShape<2>::type set_shape;


// CGP-like datatypes for coordinate sets

typedef vigra::MultiArrayShape<3>::type three_coordinate;	// coordinate in 3d MultiArray
typedef vigra::ArrayVector<three_coordinate> three_set;	// List of 3D coordinates
typedef vigra::ArrayVector<three_coordinate>::iterator three_iterator;

typedef vigra::ArrayVector<label_type> value_set;

// Datatype of parts counter in geometry file

typedef unsigned long ulong;
typedef vigra::MultiArray<1,ulong> counter_type;


#endif //CONFIGFEATURES_HXX
