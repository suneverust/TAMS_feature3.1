 /*********************************************
 * Author: Bo Sun                            *
 * Afflication: TAMS, University of Hamburg  *
 * E-Mail: bosun@informatik.uni-hamburg.de   *
 *         user_mail@QQ.com                  *
 * Date: Oct 13, 2014                        *
 * Licensing: GNU GPL license.               *
 *********************************************/

#ifndef TAMS_FEATURE_IMPL_H_
#define TAMS_FEATURE_IMPL_H_

#include "tams_feature.h"

template <typename PointInT, typename PointOuT> inline void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::setEntropyBinDim (
        size_t tams_entropy_bin_)
{
    this->tams_entropy_bin_ = tams_entropy_bin_;
}

template <typename PointInT, typename PointOuT> inline void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::setSEIAzimuthDim (
        size_t tams_sei_azimutch_dim_)
{
    this->tams_sei_azimutch_dim_ = tams_sei_azimutch_dim_;
    this->tams_sei_azimutch_spa_ = 2*M_PI/tams_sei_azimutch_dim_;
}

template <typename PointInT, typename PointOuT> inline void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::setSEIPolarDim (
        size_t tams_sei_polar_dim_)
{
    this->tams_sei_polar_dim_ = tams_sei_polar_dim_;
    this->tams_sei_polar_spa_ = M_PI/tams_sei_polar_dim_;
}

template <typename PointInT, typename PointOuT> inline int
tams::TAMSFeatureEstimation<PointInT, PointOuT>::getEntropyBinDim ()
{
    return (tams_entropy_bin_);
}

template <typename PointInT, typename PointOuT> inline int
tams::TAMSFeatureEstimation<PointInT, PointOuT>::getSEIAzimuchDim ()
{
    return (tams_sei_azimutch_dim_);
}

template <typename PointInT, typename PointOuT> inline int
tams::TAMSFeatureEstimation<PointInT, PointOuT>::getSEIPolarDim ()
{
    return (tams_sei_polar_dim_);
}

template <typename PointInT, typename PointOuT> void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::tams_cart2sph(
        float x, float y, float z,
        float& azimuth, float& polar)
{
    polar = atan2(hypot(x,y),z);
    azimuth = atan2 (y,x);
    if (azimuth<0)
        azimuth = azimuth+2*M_PI;
}

template <typename PointInT, typename PointOuT> void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::tams_vector_normalization(
        std::vector<float> &tams_vector)
{
    float max_element = (*std::max_element(tams_vector.begin(),tams_vector.end()));
    float min_element = (*std::min_element(tams_vector.begin(),tams_vector.end()));

    for (std::vector<float>::iterator itr = tams_vector.begin(); itr != tams_vector.end(); itr ++)
    {
        (*itr) = ((*itr)-min_element)/(max_element-min_element);
    }
}

template <typename PointInT, typename PointOuT> void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::tams_vector2entropy(
        const std::vector<float> & tams_vector, float &entropy)
{
    // build histogram
    std::vector<float> temp_sei_hist(tams_entropy_bin_+1);
    for (std::vector<float>::const_iterator itr = tams_vector.begin();
         itr != tams_vector.end(); itr++)
    {
        temp_sei_hist[floor((*itr)*tams_entropy_bin_)]++;
    }
    temp_sei_hist[tams_entropy_bin_-1]++;
    temp_sei_hist.pop_back();
    if(temp_sei_hist.size() != tams_entropy_bin_)
        pcl::console::print_warn(
                    "Warning: some things wrong in computing Histogram for Entropy!\n");

    //Parzen Window:[0.05,0.25,0.40,0.25,0.05]
    std::vector<float> temp_sei_hist_pad;
    temp_sei_hist_pad.push_back(0);
    temp_sei_hist_pad.push_back(0);
    for (std::vector<float>::iterator itr = temp_sei_hist.begin();
         itr != temp_sei_hist.end(); itr++)
    {
        temp_sei_hist_pad.push_back(*itr);
    }
    temp_sei_hist_pad.push_back(0);
    temp_sei_hist_pad.push_back(0);
    std::vector<float>().swap(temp_sei_hist);

    std::vector<float> tams_sei_hist;
    for (std::vector<float>::iterator itr = temp_sei_hist_pad.begin()+2;
         itr !=temp_sei_hist_pad.end()-2; itr++)
    {
        tams_sei_hist.push_back( (*(itr-2))*0.05
                                +(*(itr-1))*0.25
                                +(*itr    )*0.40
                                +(*(itr+1))*0.25
                                +(*(itr+2))*0.05);
    }
    if (tams_sei_hist.size()!=tams_entropy_bin_)
    {
        pcl::console::print_error("Error: Histogram Parzen Windows failed!\n");
        return;
    }
    std::vector<float>().swap(temp_sei_hist_pad);

    entropy = 0.0;
    for (std::vector<float>::iterator itr = tams_sei_hist.begin();
         itr!=tams_sei_hist.end(); itr++)
    {
        if ((*itr)>0)
            entropy += -(*itr)*log(*itr);
    }
    std::vector<float>().swap(tams_sei_hist);
}

template <typename PointInT, typename PointOuT> void
tams::TAMSFeatureEstimation<PointInT, PointOuT>::computeFeature(
        PointCloudOut &output)
{
    // IMPORTANT!!! OR there will be errors:
    // <1> when execute searchForNeighbors
    // <2> when user NOT set indices
    this->initCompute();

    // IMPORTANT!!! OR we could NOT assign the feature values to FeatureCloud
    output.width  = indices_->size();
    output.height = 1;
    output.points.resize(output.width*output.height);

    Eigen::Vector3f x_axis;
    Eigen::Vector3f y_axis;
    Eigen::Vector3f z_axis;

    std::vector<int> nn_indices ;
    std::vector<float> nn_dists ;

    output.is_dense = true;
    for (size_t idx = 0; idx < indices_->size(); ++idx)
    {
        if(!isFinite ((*input_).points[(*indices_)[idx]])||
                this->searchForNeighbors ((*indices_)[idx], search_parameter_, nn_indices, nn_dists)==0)
        {
            for (int d=0; d < FEATURE_SIZE; ++d)
                output.points[idx].tams_sei[d]= std::numeric_limits<float>::quiet_NaN ();

            output.is_dense = false;
            continue;
        }

        ////////////////////////////////////////////////////////////////////////////
        // The following computing LRF snippet is modified from shot_lrf.hpp
          /** \note If you use this code in any academic work, please cite:
            *
            *   - F. Tombari, S. Salti, L. Di Stefano
            *     Unique Signatures of Histograms for Local Surface Description.
            *     In Proceedings of the 11th European Conference on Computer Vision (ECCV),
            *     Heraklion, Greece, September 5-11 2010.
            *   - F. Tombari, S. Salti, L. Di Stefano
            *     A Combined Texture-Shape Descriptor For Enhanced 3D Feature Matching.
            *     In Proceedings of the 18th International Conference on Image Processing (ICIP),
            *     Brussels, Belgium, September 11-14 2011.
            *
            * \author Samuele Salti, Federico Tombari
            * \ingroup features
            */
        const Eigen::Vector4f& central_point = (*input_).points[(*indices_)[idx]].getVector4fMap();

        Eigen::Matrix<double, Eigen::Dynamic, 4> vij (nn_indices.size(), 4);
        Eigen::Matrix3d cov_m = Eigen::Matrix3d::Zero();
        double distance = 0.0;
        double sum = 0.0;
        int valid_nn_points = 0;

        for (size_t i_idx = 0; i_idx < nn_indices.size(); ++i_idx)
        {
            Eigen::Vector4f pt = (*surface_).points[nn_indices[i_idx]].getVector4fMap();
            if(pt.head<3>()==central_point.head<3>())
                continue;
            // Difference between current point and central_point.head<3>()
            vij.row(valid_nn_points).matrix() = (pt-central_point).cast<double>();
            vij(valid_nn_points,3) = 0;

            distance = search_parameter_ - sqrt(nn_dists[i_idx]);

            // Multiply vij*vij'
            cov_m += distance * (vij.row(valid_nn_points).head<3>().transpose()*
                                 vij.row(valid_nn_points).head<3>());

            sum += distance;
            valid_nn_points++;
        }
        if (valid_nn_points < 5)
        {
            for (int d=0; d < FEATURE_SIZE; ++d)
                output.points[idx].tams_sei[d]= std::numeric_limits<float>::quiet_NaN ();

            output.is_dense = false;
            continue;
        }

        cov_m /=sum;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver (cov_m);
        const double& e1c = solver.eigenvalues()[0];
        const double& e2c = solver.eigenvalues()[1];
        const double& e3c = solver.eigenvalues()[2];

        if (!pcl_isfinite(e1c) || !pcl_isfinite(e2c) || !pcl_isfinite(e3c))
        {
            for (int d=0; d < FEATURE_SIZE; ++d)
                output.points[idx].tams_sei[d]= std::numeric_limits<float>::quiet_NaN ();

            output.is_dense = false;
            continue;
        }

        // Disambiguation
        Eigen::Vector4d v1 = Eigen::Vector4d::Zero();
        Eigen::Vector4d v3 = Eigen::Vector4d::Zero();
        v1.head<3>().matrix() = solver.eigenvectors().col(2);
        v3.head<3>().matrix() = solver.eigenvectors().col(0);

        int plusNormal = 0, plusTangentDirection1 = 0;
        for (int ne = 0; ne < valid_nn_points; ne++)
        {
            double dp = vij.row(ne).dot(v1);
            if (dp >= 0)
                plusTangentDirection1++;

            dp = vij.row(ne).dot(v3);
            if (dp >=0)
                plusNormal++;
        }

        //TANGENT
        plusTangentDirection1 = 2*plusTangentDirection1 - valid_nn_points;
        if (plusTangentDirection1 == 0)
        {
              int points = 5;
              int medianIndex = valid_nn_points/2;

              for (int i = -points/2; i <= points/2; i++)
                  if ( vij.row (medianIndex - i).dot (v1) > 0)
                      plusTangentDirection1 ++;

              if (plusTangentDirection1 < points/2+1)
                  v1 *= - 1;
          }
        else if (plusTangentDirection1 < 0)
          v1 *= - 1;

        //Normal
        plusNormal = 2*plusNormal - valid_nn_points;
        if (plusNormal == 0)
        {
              int points = 5;
              int medianIndex = valid_nn_points/2;

              for (int i = -points/2; i <= points/2; i++)
                  if ( vij.row (medianIndex - i).dot (v3) > 0)
                      plusNormal ++;

              if (plusNormal < points/2+1)
                  v3 *= - 1;
          }
        else if (plusNormal < 0)
          v3 *= - 1;

        x_axis = v1.head<3>().cast<float>();
        z_axis = v3.head<3>().cast<float>();
        y_axis = z_axis.cross(x_axis);

        ////////////////////////////////////////////////////////////////////////////

        // Check the local reference frame
        if(fabs(x_axis[0]*y_axis[0]+x_axis[1]*y_axis[1]+x_axis[2]*y_axis[2])>1E-4f)
        {
            pcl::console::print_warn("The LRF of point: idx = %d got some warning\n", idx);
            pcl::console::print_warn("x.cross(y) is %f",
                                     x_axis[0]*y_axis[0]+x_axis[1]*y_axis[1]+x_axis[2]*y_axis[2]);
        }

        if(fabs(x_axis[0]*z_axis[0]+x_axis[1]*z_axis[1]+x_axis[2]*z_axis[2])>1E-4f)
        {
            pcl::console::print_warn("The LRF of point: idx = %d got some warning\n", idx);
            pcl::console::print_warn("x.cross(z) is %f",
                                     x_axis[0]*z_axis[0]+x_axis[1]*z_axis[1]+x_axis[2]*z_axis[2]);
        }
        if(fabs(y_axis[0]*z_axis[0]+y_axis[1]*z_axis[1]+y_axis[2]*z_axis[2])>1E-4f)
        {
            pcl::console::print_warn("The LRF of point: idx = %d got some warning\n", idx);
            pcl::console::print_warn("y.cross(z) is %f",
                                     y_axis[0]*z_axis[0]+y_axis[1]*z_axis[1]+y_axis[2]*z_axis[2]);
        }

        // Main implementation
        // Point Division
        float temp_az, temp_polar;
        size_t temp_sei_azth, temp_sei_polarth;

        Eigen::Array<std::vector<float>, Eigen::Dynamic, Eigen::Dynamic>
                TAMS_sei_points(tams_sei_azimutch_dim_, tams_sei_polar_dim_);
        for (size_t Neighbor_idx = 0; Neighbor_idx < nn_indices.size(); Neighbor_idx++)
        {
            if (pcl::utils::equal(nn_dists[Neighbor_idx], 0.0f, 1E-4f))
                continue;
            Eigen::Vector3f neighbor = (*surface_).points[nn_indices[Neighbor_idx]].getVector3fMap();
            // compute neighbour polar coordinated in current LRF
            Eigen::Vector3f proj;
            pcl::geometry::project(neighbor, central_point.head<3>(), z_axis, proj);
            proj = proj - central_point.head<3>();
            proj.normalize();
            // compute the angle between the projection and the x_axis in [0, 2*M_PI]
            if (!pcl_isfinite(proj(0))|!pcl_isfinite(proj(1))|!pcl_isfinite(proj(2)))
                temp_az = 0;
            else
            {
                Eigen::Vector3f cross = x_axis.cross(proj);
                temp_az = std::atan2(cross.norm(), x_axis.dot(proj));
                temp_az = cross.dot(z_axis) < 0.0f ? (2*M_PI-temp_az):temp_az;
            }
            // compute the angle between the neighbor and the z_axis (normal) in [0, M_PI]
            Eigen::Vector3f neig = neighbor - central_point.head<3>();
            neig.normalize();
            temp_polar = z_axis.dot(neig);
            temp_polar = acosf(std::min(1.0f, std::max(-1.0f, temp_polar)));

            if (temp_az < tams_sei_azimutch_spa_/2 ||
                    temp_az >= 2*M_PI-tams_sei_azimutch_spa_/2)
                temp_sei_azth = 0;
            else
                temp_sei_azth = floor((temp_az-tams_sei_azimutch_spa_/2)/tams_sei_azimutch_spa_)+1;
            if (temp_polar-M_PI < 1E-4f)
                temp_sei_polarth = tams_sei_polar_dim_-1;
            else
                temp_sei_polarth = floor(temp_polar/tams_sei_polar_spa_);

            if (temp_sei_azth >= tams_sei_azimutch_dim_)
                std::cout << "temp_az"<<temp_az << std::endl;
            if (temp_sei_polarth >= tams_sei_polar_dim_)
                std::cout << "temp_polar "<< temp_polar << std::endl;
            TAMS_sei_points (temp_sei_azth, temp_sei_polarth).push_back (nn_dists[Neighbor_idx]);
        }

        // Entropy Computation
        std::vector<float> temp_sei_points;
        Eigen::MatrixXf TAMS_sei = Eigen::MatrixXf::Zero(tams_sei_azimutch_dim_, tams_sei_polar_dim_);
        for (temp_sei_azth = 0; temp_sei_azth < tams_sei_azimutch_dim_; temp_sei_azth++)
        {
            for (temp_sei_polarth = 0; temp_sei_polarth < tams_sei_polar_dim_; temp_sei_polarth++)
            {
                temp_sei_points = TAMS_sei_points (temp_sei_azth, temp_sei_polarth);

                if (temp_sei_points.size() < 5)
                {
                    TAMS_sei(temp_sei_azth, temp_sei_polarth) = 0;
                    continue;
                }

                if ((*std::max_element(temp_sei_points.begin(),temp_sei_points.end())) ==
                        (*std::min_element(temp_sei_points.begin(),temp_sei_points.end())))
                {
                    TAMS_sei(temp_sei_azth, temp_sei_polarth) = 0;
                    continue;
                }

                this->tams_vector_normalization(temp_sei_points);

                this->tams_vector2entropy(temp_sei_points,TAMS_sei(temp_sei_azth, temp_sei_polarth));
            }
        }
        std::vector<float>().swap(temp_sei_points);

        // Entropy Resize
        Eigen::VectorXf TAMS_sei_real = Eigen::VectorXf::Zero(
                    tams_sei_azimutch_dim_*tams_sei_polar_dim_);
        TAMS_sei.resize(tams_sei_azimutch_dim_*tams_sei_polar_dim_,1);
        TAMS_sei_real << TAMS_sei;
        TAMS_sei.resize(0,0);

        for (int i = 0; i < TAMS_sei_real.size();i++)
            output.points[idx].tams_sei[i]=TAMS_sei_real(i);
    }
}
#endif/*TAMS_FEATURE_IMPL_H_*/

