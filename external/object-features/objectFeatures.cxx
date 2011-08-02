#include "objectFeatures.hxx"

namespace features
{

    feature_array extractVolume(three_set &coordinates, value_set &intensities)
    {
        feature_array ret (array_shape(1),double(coordinates.size()));
        return ret;
    }


    feature_array extractPosition(three_set &coordinates, value_set &intensities)
    {
        int size = coordinates.size();
        feature_array ret (array_shape(12),-1.);

        //catch cases with empty or small lists
        if(size == 0){
            return ret;
        }
        if(size == 1){
            ret[0] = coordinates[0][0];
            ret[1] = coordinates[0][1];
            ret[2] = coordinates[0][2];
            return ret;
        }

        long int sum_x = 0, sum_y = 0, sum_z = 0;
        double var_x = 0, var_y = 0, var_z = 0;
        double skew_x = 0, skew_y = 0, skew_z = 0;
        double kurt_x = 0, kurt_y = 0, kurt_z = 0;

        // calculate mean first
        for(three_iterator it = coordinates.begin(); it != coordinates.end(); ++it){
            sum_x += (*it)[0];
            sum_y += (*it)[1];
            sum_z += (*it)[2];
        }

        ret[0] = double(sum_x)/size;
        ret[1] = double(sum_y)/size;
        ret[2] = double(sum_z)/size;

        // then calculate central moments
        for(three_iterator it = coordinates.begin(); it != coordinates.end(); ++it){
            double delta_x = (*it)[0] - ret[0];
            double delta_y = (*it)[1] - ret[1];
            double delta_z = (*it)[2] - ret[2];
            var_x += std::pow(delta_x,2);
            var_y += std::pow(delta_y,2);
            var_z += std::pow(delta_z,2);

            //*
            skew_x += std::pow(delta_x,3);
            skew_y += std::pow(delta_y,3);
            skew_z += std::pow(delta_z,3);

            kurt_x += std::pow(delta_x,4);
            kurt_y += std::pow(delta_y,4);
            kurt_z += std::pow(delta_z,4);
            //*/
        }

        ret[3] = var_x/double(size);
        ret[4] = var_y/double(size);
        ret[5] = var_z/double(size);

        //*
        if(ret[3] != 0){
            ret[6] = skew_x/double(size)/std::pow(double(ret[3]),3./2.);
            ret[9] = kurt_x/double(size)/std::pow(ret[3],2) - 3;
        }
        if(ret[4] != 0){
            ret[7] = skew_y/double(size)/std::pow(double(ret[4]),3./2.);
            ret[10] = kurt_y/double(size)/std::pow(ret[4],2) - 3;
        }
        if(ret[5] != 0){
            ret[8] = skew_z/double(size)/std::pow(double(ret[5]),3./2.);
            ret[11] = kurt_z/double(size)/std::pow(ret[5],2) - 3;
        }
        //*/
        return ret;
    }


    feature_array extractWeightedPosition(three_set &coordinates, value_set &intensities)
    {
        int size = coordinates.size();

        long int sum_x = 0, sum_y = 0, sum_z = 0;
        double var_x = 0, var_y = 0, var_z = 0;
        double skew_x = 0, skew_y = 0, skew_z = 0;
        double kurt_x = 0, kurt_y = 0, kurt_z = 0;

        int sum_intensity = 0;

        feature_array ret (array_shape(12),-1.);
        if(size == 0){
            return ret;
        }

        // calculate mean first
        for(int i = 0; i < size; i++){
            int intensity = intensities[i];
            sum_x += intensity*coordinates[i][0];
            sum_y += intensity*coordinates[i][1];
            sum_z += intensity*coordinates[i][2];
            sum_intensity += intensity;
        }

        if(sum_intensity == 0){
            return ret;
        }

        ret[0] = double(sum_x)/double(sum_intensity);
        ret[1] = double(sum_y)/double(sum_intensity);
        ret[2] = double(sum_z)/double(sum_intensity);

        if(size <= 1){
            return ret;
        }

        // then calculate central moments
        for(int i = 0; i < size; i++){
            int intensity = intensities[i];
            double delta_x = coordinates[i][0] - ret[0];
            double delta_y = coordinates[i][1] - ret[1];
            double delta_z = coordinates[i][2] - ret[2];
            var_x += intensity*std::pow(delta_x,2);
            var_y += intensity*std::pow(delta_y,2);
            var_z += intensity*std::pow(delta_z,2);

            //*
            skew_x += intensity*std::pow(delta_x,3);
            skew_y += intensity*std::pow(delta_y,3);
            skew_z += intensity*std::pow(delta_z,3);

            kurt_x += intensity*std::pow(delta_x,4);
            kurt_y += intensity*std::pow(delta_y,4);
            kurt_z += intensity*std::pow(delta_z,4);
            //*/
        }

        ret[3] = var_x/double(sum_intensity);
        ret[4] = var_y/double(sum_intensity);
        ret[5] = var_z/double(sum_intensity);

        //*
        if(ret[3] != 0){
            ret[6] = skew_x/double(sum_intensity)/std::pow(double(ret[3]),3./2.);
            ret[9] = kurt_x/double(sum_intensity)/std::pow(ret[3],2) - 3;
        }
        if(ret[4] != 0){
            ret[7] = skew_y/double(sum_intensity)/std::pow(double(ret[4]),3./2.);
            ret[10] = kurt_y/double(sum_intensity)/std::pow(ret[4],2) - 3;
        }
        if(ret[5] != 0){
            ret[8] = skew_z/double(sum_intensity)/std::pow(double(ret[5]),3./2.);
            ret[11] = kurt_z/double(sum_intensity)/std::pow(ret[5],2) - 3;
        }
        //*/

        return ret;
    }


    feature_array extractPrincipalComponents(three_set &coordinates, value_set &intensities)
    {
        int size = coordinates.size();
        feature_array ret (array_shape(12),-1);

        if(size <= 1){
            return ret;
        }

        // set up vigra matrix with coordinates shifted to zero mean
        feature_array mv = extractPosition(coordinates, intensities);
        vigra::linalg::Matrix<double> coord_matrix (size,3,0.);
        for(int i = 0; i < size; i++){
            coord_matrix(i,0) = double(coordinates[i][0] - mv[0]);
            coord_matrix(i,1) = double(coordinates[i][1] - mv[1]);
            coord_matrix(i,2) = double(coordinates[i][2] - mv[2]);
        }

        // claculate covariance matrix
        vigra::linalg::Matrix<double> cov_matrix (3,3,0.);
        vigra::linalg::covarianceMatrixOfColumns(coord_matrix, cov_matrix);

        // calculate eigenvalues and eigenvectors
        vigra::linalg::Matrix<double> ev_matrix (3,3,0.);	//matrix for eigenvectors
        vigra::linalg::Matrix<double> ew_vector (3,1,0.);	//matrix for eigenvalues

        vigra::linalg::symmetricEigensystem(cov_matrix,ew_vector,ev_matrix);

        // first fill the eigenvalues to the output vector
        ret[0] = ew_vector(0,0);
        ret[1] = ew_vector(1,0);
        ret[2] = ew_vector(2,0);

        // then the eigenvectors
        ret[3]  = ev_matrix(0,0);
        ret[4]  = ev_matrix(1,0);
        ret[5]  = ev_matrix(2,0);
        ret[6]  = ev_matrix(0,1);
        ret[7]  = ev_matrix(1,1);
        ret[8]  = ev_matrix(2,1);
        ret[9]  = ev_matrix(0,2);
        ret[10] = ev_matrix(1,2);
        ret[11] = ev_matrix(2,2);

        return ret;
    }


    feature_array extractBoundingBox(three_set &coordinates, value_set &intensities)
    {
        int size = coordinates.size();
        feature_array ret (array_shape(7),-1);

        //catch cases with empty or small lists
        if(size == 0){
            return ret;
        }

        // init with some real data value
        ret[0] = coordinates[0][0];	// --> x Min
        ret[3] = coordinates[0][0];	// --> x Max
        ret[1] = coordinates[0][1];	// --> y Min
        ret[4] = coordinates[0][1];	// --> y Max
        ret[2] = coordinates[0][2];	// --> z Min
        ret[5] = coordinates[0][2];	// --> z Max

        // go through list and find extremal cordinates
        for(three_iterator it = coordinates.begin(); it != coordinates.end(); ++it){
            if( (*it)[0] < ret[0] ){
                ret[0] = (*it)[0];
            }
            if( (*it)[0] > ret[3] ){
                ret[3] = (*it)[0];
            }
            if( (*it)[1] < ret[1] ){
                ret[1] = (*it)[1];
            }
            if( (*it)[1] > ret[4] ){
                ret[4] = (*it)[1];
            }
            if( (*it)[2] < ret[2] ){
                ret[2] = (*it)[2];
            }
            if( (*it)[2] > ret[5] ){
                ret[5] = (*it)[2];
            }
        }

        // sizes along dimensions
        double sx = ret[3] - ret[0] + 1;
        double sy = ret[4] - ret[1] + 1;
        double sz = ret[5] - ret[2] + 1;

        // calculate fill factor
        if(sx!=0 && sy!=0 && sz!=0){
            ret[6] = double(size) / (sx * sy * sz);
        }

        return ret;
    }



    feature_array extractIntensity(three_set &coordinates, value_set &intensities)
    {
        feature_array ret (array_shape(4),-1);
        int size = intensities.size();

        //catch cases with empty or small lists
        if(size == 0){
            return ret;
        }
        if(size == 1){
            ret[0] = intensities[0];
            return ret;
        }

        //calculate mean first

        double sum_intensity = 0.;

        for(int i = 0; i < size; i++){
            sum_intensity += double(intensities[i]);
        }

        ret[0] = sum_intensity / double(size);


        // Calculate variance, skew and kurtosis
        double var = 0;
        double skew = 0;
        double kurt = 0;

        for(int i = 0; i < size; i++){
            double delta = double(intensities[i]) - ret[0];
            var += std::pow(delta,2);
            skew += std::pow(delta,3);
            kurt += std::pow(delta,4);
        }

        ret[1] = var / double(size);

        if(ret[1] != 0){
            ret[2] = skew/double(size)/std::pow(double(ret[1]),3./2.);
            ret[3] = kurt/double(size)/std::pow(ret[1],2) - 3;
        }

        return ret;
    }



    feature_array extractMinMaxIntensity(three_set &coordinates, value_set &intensities)
    {
        feature_array ret (array_shape(9),-1);
        int size = intensities.size();

        //catch cases with empty or small lists
        if(size == 0){
            return ret;
        }

        std::vector<feature_type> intensityList;
        for(int i = 0; i < size; i++)
        {
            intensityList.push_back(intensities[i]);
        }

        std::sort(intensityList.begin(),intensityList.end());

        // minimum/maximum intensity
        ret[0] = intensityList.front();
        ret[1] = intensityList.back();

        // find the 5%, 10%, 20%, 50%, 80%, 90%, 95% quanitles
        int i = 0;

        // 5% quantile
        while(i < 0.05*size && i < size-1){
            i++;
        }

        ret[2] = intensityList[i];

        // 10% quantile
        while(i < 0.1*size && i < size-1){
            i++;
        }

        ret[3] = intensityList[i];

        // 20% quantile
        while(i < 0.2*size && i < size-1){
            i++;
        }

        ret[4] = intensityList[i];

        // 50% quantile
        while(i < 0.5*size && i < size-1){
            i++;
        }

        ret[5] = intensityList[i];

        // 80% quantile
        while(i < 0.8*size && i < size-1){
            i++;
        }

        ret[6] = intensityList[i];

        // 90% quantile
        while(i < 0.9*size && i < size-1){
            i++;
        }

        ret[7] = intensityList[i];

        // 95% quantile
        while(i < 0.95*size && i < size-1){
            i++;
        }

        ret[8] = intensityList[i];

        return ret;
    }



    feature_array extractMaxIntensity(three_set &coordinates, value_set &intensities)
    {
        feature_array ret (array_shape(4),-1);
        int size = intensities.size();

        feature_array com = extractWeightedPosition(coordinates,intensities);

        //catch cases with empty or small lists
        if(size == 0){
            return ret;
        }

        int max = intensities[0];
        int x1 = coordinates[0][0];
        int x2 = coordinates[0][1];
        int x3 = coordinates[0][2];

        for(int i = 1; i < size; i++)
        {
            if(intensities[i] > max){

                max = intensities[i];
                x1 = coordinates[i][0];
                x2 = coordinates[i][1];
                x3 = coordinates[i][2];

            }else if(intensities[i] == max){

                double d_old = std::sqrt(std::pow(x1-com[0],2)+std::pow(x2-com[1],2)+std::pow(x3-com[2],2));
                double d_new = std::sqrt(std::pow(coordinates[i][0]-com[0],2)+std::pow(coordinates[i][1]-com[1],2)+std::pow(coordinates[i][2]-com[2],2));
                if(d_new < d_old){
                    max = intensities[i];
                    x1 = coordinates[i][0];
                    x2 = coordinates[i][1];
                    x3 = coordinates[i][2];
                }

            }

        }

        ret[0] = max;
        ret[1] = x1;
        ret[2] = x2;
        ret[3] = x3;

        return ret;
    }



    feature_array extractPairwise(three_set &coordinates, value_set &intensities)
    {
        feature_array ret (array_shape(4),0.);
        int size = intensities.size();

        //catch cases with empty or small lists
        if(size <= 1){
            return ret;
        }

        //fill a zero-padded volume
        feature_array bbox (array_shape(7),0.);
        bbox = extractBoundingBox(coordinates,intensities);

        vigra::MultiArray<3,feature_type> volume (volume_shape(bbox[3]-bbox[0]+1,bbox[4]-bbox[1]+1,bbox[5]-bbox[2]+1),0.);

        for(int i = 0; i < size; i++)
        {
            coordinate_type x0 = coordinates[i][0] - bbox[0];
            coordinate_type x1 = coordinates[i][1] - bbox[1];
            coordinate_type x2 = coordinates[i][2] - bbox[2];
            volume(x0,x1,x2) = intensities[i];
        }

        //go over volume and sum up pairwise costs
        for(int i = 0; i < volume.size(0); i++){
            for(int j = 0; j < volume.size(1); j++){
                for(int k = 0; k < volume.size(2); k++){

                    if(i+1 < volume.size(0)){
                        // relative sum over absolute distances
                        ret[0] += std::abs(volume(i+1,j,k) - volume(i,j,k))/double(size);
                        // relative sum over squared distances
                        ret[1] += std::pow(volume(i+1,j,k) - volume(i,j,k),2)/double(size);
                        if(i>0){
                            // sum over symmetric derivatives
                            ret[2] += std::abs(volume(i+1,j,k) - volume(i-1,j,k))/double(size);
                            // second derivative
                            ret[3] += std::abs(2*volume(i,j,k) - volume(i-1,j,k) - volume(i+1,j,k))/double(size);
                        }
                    }
                    if(j+1 < volume.size(1)){
                        ret[0] += std::abs(volume(i,j+1,k) - volume(i,j,k))/double(size);
                        ret[1] += std::pow(volume(i,j+1,k) - volume(i,j,k),2)/double(size);
                        if(j>0){
                            ret[2] += std::abs(volume(i,j+1,k) - volume(i,j-1,k))/double(size);
                            ret[3] += std::abs(2*volume(i,j,k) - volume(i,j-1,k) - volume(i,j+1,k))/double(size);
                        }
                    }
                    if(k+1 < volume.size(2)){
                        ret[0] += std::abs(volume(i,j,k+1) - volume(i,j,k))/double(size);
                        ret[1] += std::pow(volume(i,j,k+1) - volume(i,j,k),2)/double(size);
                        if(k>0){
                            ret[2] += std::abs(volume(i,j,k+1) - volume(i,j,k-1))/double(size);
                            ret[3] += std::abs(2*volume(i,j,k) - volume(i,j,k-1) - volume(i,j,k+1))/double(size);
                        }
                    }
                }
            }
        }

        return ret;
    }



    feature_array extractHOG(three_set &coordinates, value_set &intensities)
    {
        feature_array ret (array_shape(20),0.);
        int size = intensities.size();

        //catch cases with empty or small lists
        if(size <= 1){
            return ret;
        }

        //fill a zero-padded volume
        feature_array bbox (array_shape(7),0.);
        bbox = extractBoundingBox(coordinates,intensities);

        vigra::MultiArray<3,feature_type> volume (volume_shape(bbox[3]-bbox[0]+1,bbox[4]-bbox[1]+1,bbox[5]-bbox[2]+1),0.);

        for(int i = 0; i < size; i++)
        {
            coordinate_type x0 = coordinates[i][0] - bbox[0];
            coordinate_type x1 = coordinates[i][1] - bbox[1];
            coordinate_type x2 = coordinates[i][2] - bbox[2];
            volume(x0,x1,x2) = intensities[i];
        }

        //set up the 20 binning directions
        double phi = (1. + std::sqrt(5.))/2.;
        double threshold = 1.29107;
        vigra::MultiArray<2,feature_type> icodir (matrix_shape(20,3));
        icodir(0,0) =  1; icodir(0,1) =  1; icodir(0,2) =  1;
        icodir(1,0) =  1; icodir(1,1) =  1; icodir(1,2) = -1;
        icodir(2,0) =  1; icodir(2,1) = -1; icodir(2,2) =  1;
        icodir(3,0) =  1; icodir(3,1) = -1; icodir(3,2) = -1;
        icodir(4,0) = -1; icodir(4,1) =  1; icodir(4,2) =  1;
        icodir(5,0) = -1; icodir(5,1) =  1; icodir(5,2) = -1;
        icodir(6,0) = -1; icodir(6,1) = -1; icodir(6,2) =  1;
        icodir(7,0) = -1; icodir(7,1) = -1; icodir(7,2) = -1;

        icodir(8,0) =  0; icodir(8,1) =   1/phi; icodir(8,2) =   phi;
        icodir(9,0) =  0; icodir(9,1) =   1/phi; icodir(9,2) =  -phi;
        icodir(10,0) = 0; icodir(10,1) = -1/phi; icodir(10,2) =  phi;
        icodir(11,0) = 0; icodir(11,1) = -1/phi; icodir(11,2) = -phi;

        icodir(12,0) =  1/phi; icodir(12,1) =  phi; icodir(12,2) = 0;
        icodir(13,0) =  1/phi; icodir(13,1) = -phi; icodir(13,2) = 0;
        icodir(14,0) = -1/phi; icodir(14,1) =  phi; icodir(14,2) = 0;
        icodir(15,0) = -1/phi; icodir(15,1) = -phi; icodir(15,2) = 0;

        icodir(16,0) =  phi; icodir(16,1) = 0; icodir(16,2) =  1/phi;
        icodir(17,0) =  phi; icodir(17,1) = 0; icodir(17,2) = -1/phi;
        icodir(18,0) = -phi; icodir(18,1) = 0; icodir(18,2) =  1/phi;
        icodir(19,0) = -phi; icodir(19,1) = 0; icodir(19,2) = -1/phi;

        //go over volume and collect votes.
        for(int i = 1; i < volume.size(0)-1; i++){
            for(int j = 1; j < volume.size(1)-1; j++){
                for(int k = 1; k < volume.size(2)-1; k++){
                    feature_type grad_x = volume(i+1,j,k) - volume(i-1,j,k);
                    feature_type grad_y = volume(i,j+1,k) - volume(i,j-1,k);
                    feature_type grad_z = volume(i,j,k+1) - volume(i,j,k-1);
                    feature_type norm = std::sqrt(std::pow(grad_x,2)+std::pow(grad_y,2)+std::pow(grad_z,2));
                    if(norm > 0){
                        vigra::ArrayVector<feature_type> votes (20, 0.);
                        for(int s = 0; s < 20; s++){
                            feature_type v = (grad_x*icodir(s,0) + grad_y*icodir(s,1) + grad_z*icodir(s,2))/norm;
                            if(v > threshold){
                                votes[s] = v;
                            }
                        }
                        double sum_votes = 0.;
                        for(int s = 0; s < 20; s++){
                            sum_votes += votes[s];
                        }
                        if(sum_votes <= 0)
                            std::cout << "error: sum <= 0" << std::endl;
                        for(int s = 0; s < 20; s++){
                            ret[s] += norm*votes[s]/sum_votes;
                        }
                    }
                }
            }
        }

        for(int s = 0; s < 20; s++){
            ret[s] /= double(size);
        }

        return ret;
    }



    feature_array extractSGF(three_set &coordinates, value_set &intensities)
    {

        feature_array ret (array_shape(48),0.);
        int size = intensities.size();

        //catch cases with empty or small lists
        if(size <= 1){
            return ret;
        }

        //get minimum/maximum coordinates
        feature_array bbox (array_shape(7),0.);
        bbox = extractBoundingBox(coordinates,intensities);

        // fill a volume with intensity values
        vigra::MultiArray<3,raw_type> volume (volume_shape(bbox[3]-bbox[0]+1,bbox[4]-bbox[1]+1,bbox[5]-bbox[2]+1),0.);

        for(int i = 0; i < size; i++)
        {
            coordinate_type x0 = coordinates[i][0] - bbox[0];
            coordinate_type x1 = coordinates[i][1] - bbox[1];
            coordinate_type x2 = coordinates[i][2] - bbox[2];
            volume(x0,x1,x2) = raw_type(intensities[i]);
        }

        // Calculate the Features
        return SGF::SGFeatures<3,raw_type>(volume, 24);

        return ret;
    }



} /* namespace features */

