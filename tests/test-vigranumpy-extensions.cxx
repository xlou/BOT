#include "objectFeatures.hxx"
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    three_set coordinates;
    value_set intensities;

    // create an object
    for(int x = 10; x < 15; x++) {
        for(int y = 20; y < 30; y++) {
            coordinates.push_back( three_coordinate(x, y, 40));
            if (x >= 12 && x <= 13 && y >= 23 && y <= 26)
                intensities.push_back(2);
            else
                intensities.push_back(1);
        }
    }

    // volume
    {
        features::ObjectVolume<unsigned short, 3 > o(coordinates, intensities);
        feature_array volume_obj = o.get();
        std::cout << "volume: " << 
            volume_obj[0] << std::endl;
    }

    // position
    {
        features::ObjectPosition<unsigned short,3 > o(coordinates,intensities);
        feature_array position_obj = o.get();
        std::cout << "position: " << 
            position_obj[0] << ", " << position_obj[1] << ", " << position_obj[2] << std::endl;
    }

    // bounding box
    {
        features::ObjectBoundingBox<unsigned short, 3 > o(coordinates, intensities);
        feature_array bbox_obj = o.get();
        std::cout << "bounding box: " << 
            bbox_obj[0] << ", " << bbox_obj[1] << ", " << bbox_obj[2] << " - " << 
            bbox_obj[3] << ", " << bbox_obj[4] << ", " << bbox_obj[5] << std::endl;
    }

    // principal component
    {
        features::ObjectPrincipalComponents<unsigned short, 3 > o(coordinates, intensities);
        feature_array pc_obj = o.get();
        std::cout << "principal components: " << 
            pc_obj[0] << ", " << pc_obj[1] << ", " << pc_obj[2] << std::endl;
    }

    // eccentricity
    {
        features::ObjectPrincipalComponents<unsigned short, 3 > o(coordinates, intensities);
        feature_array pc_obj = o.get();
        std::cout << "eccentricity: " << 
            pc_obj[0]/pc_obj[1] << std::endl;
    }

    // intensity mean/deviation
    {
        features::ObjectIntensity<unsigned short,3 > o(coordinates, intensities);
        feature_array fint_obj = o.get();
        std::cout << "intensity mean/deviation/sum: " << 
            fint_obj[0] << ", " << sqrt(fint_obj[1]) << ", " << fint_obj[4] << std::endl;
    }

    return 0;
}
