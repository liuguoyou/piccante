 /*

PICCANTE
The hottest HDR imaging library!
http://piccantelib.net

Copyright (C) 2014
Visual Computing Laboratory - ISTI CNR
http://vcg.isti.cnr.it
First author: Francesco Banterle

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

*/

#include <QCoreApplication>

//This means that OpenGL acceleration layer is disabled
#define PIC_DISABLE_OPENGL

#define EIGEN_DONT_VECTORIZE

#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

#include "piccante.hpp"

Eigen::Vector2i Projection(Eigen::Matrix34d &M, Eigen::Vector3d p, double cx, double cy, double fx, double fy, double lambda)
{
    Eigen::Vector4d p_t = Eigen::Vector4d(p[0], p[1], p[2], 1.0);
    Eigen::Vector2i out;
    Eigen::Vector3d proj = M * p_t;
    proj[0] /= proj[2];
    proj[1] /= proj[2];

    double x_cx =  (proj[0] - cx);
    double y_cy =  (proj[1] - cy);

    double dx = x_cx / fx;
    double dy = y_cy / fy;
    double rho_sq = dx * dx + dy * dy;

    double factor = 1.0 / (1.0 + rho_sq * lambda);

    proj[0] = x_cx * factor + cx;
    proj[1] = y_cy * factor + cy;

    out[0] = int(proj[0]);
    out[1] = int(proj[1]);

    return out;
}

int main(int argc, char *argv[])
{
    Q_UNUSED(argc);
    Q_UNUSED(argv);

    printf("Reading an LDR images...");

    //estimating K matrix from camera
    double fx = pic::getFocalLengthPixels(18.0, 22.3, 2592.0);
    double fy = pic::getFocalLengthPixels(18.0, 14.9, 1728.0);
    Eigen::Matrix3d K = pic::getIntrinsicsMatrix(fx, fy, 2592.0 / 2.0, 1728.0 / 2.0);

    pic::Image img0, img1;
    img0.Read("../data/input/triangulation/venice_campo_s_polo_left.jpg", pic::LT_NOR);
    img1.Read("../data/input/triangulation/venice_campo_s_polo_right.jpg", pic::LT_NOR);

    printf("Ok\n");

    printf("Are they both valid?");
    if(img0.isValid() && img1.isValid()) {
        printf("OK\n");

        //output corners
        std::vector< Eigen::Vector3f > corners_from_img0;
        std::vector< Eigen::Vector3f > corners_from_img1;

        //computing the luminance images
        pic::Image *L0 = pic::FilterLuminance::Execute(&img0, NULL, pic::LT_CIE_LUMINANCE);
        pic::Image *L1 = pic::FilterLuminance::Execute(&img1, NULL, pic::LT_CIE_LUMINANCE);

        //getting corners
        printf("Extracting corners...\n");
        pic::HarrisCornerDetector hcd(2.5f, 5);
        hcd.Compute(L0, &corners_from_img0);
        hcd.Compute(L1, &corners_from_img1);

        //computing ORB descriptors for each corner and image
        //Computing luminance images
        pic::Image *L0_flt = pic::FilterGaussian2D::Execute(L0, NULL, 2.5f);
        pic::Image *L1_flt = pic::FilterGaussian2D::Execute(L1, NULL, 2.5f);

        printf("Computing ORB descriptors...\n");

        pic::ORBDescriptor b_desc(31, 512);

        std::vector< unsigned int *> descs0;
        for(unsigned int i=0; i<corners_from_img0.size(); i++) {
            int x = corners_from_img0[i][0];
            int y = corners_from_img0[i][1];

            descs0.push_back(b_desc.get(L0_flt, x, y));
        }

        std::vector< unsigned int *> descs1;
        for(unsigned int i=0; i<corners_from_img1.size(); i++) {
            int x = corners_from_img1[i][0];
            int y = corners_from_img1[i][1];

            descs1.push_back(b_desc.get(L1_flt, x, y));
        }

        printf("Matching ORB descriptors...\n");
        std::vector< Eigen::Vector3i > matches;

        int n = b_desc.getDescriptorSize();

        printf("Descriptor size: %d\n", n);

        //pic::BinaryFeatureBruteForceMatcher bffm_bin(&descs1, n);
        pic::BinaryFeatureLSHMatcher bffm_bin(&descs1, n);

        bffm_bin.getAllMatches(&descs0, matches);

        printf("Matches:\n");
        std::vector< Eigen::Vector2f > m0, m1;

        for(unsigned int i=0; i<matches.size(); i++) {
            int I0 = matches[i][0];
            int I1 = matches[i][1];

            Eigen::Vector2f x, y;
            x[0] = corners_from_img0[I0][0];
            x[1] = corners_from_img0[I0][1];

            y[0] = corners_from_img1[I1][0];
            y[1] = corners_from_img1[I1][1];

            m0.push_back(x);
            m1.push_back(y);

            printf("I1: %d (%d %d) -- I2: %d (%d %d) -- Score: %d\n", I0, int(x[0]), int(x[1]), I1, int(y[0]), int(y[1]), matches[i][2]);
        }

        printf("\n Total matches: (%lu | %lu)\n", m0.size(), m1.size());

        printf("\nEstimating the fundamental matrix F from the matches...");
        std::vector< unsigned int > inliers;
        Eigen::Matrix3d F = pic::estimateFundamentalRansac(m0, m1, inliers, 1000000, 0.5);

        //non-linear refinement using Nelder-Mead        
        pic::NelderMeadOptFundamental nmf(m0, m1, inliers);

        float F_data_opt[9];
        nmf.run(pic::getLinearArrayFromMatrix(F), 9, 1e-4f, 10000, &F_data_opt[0]);
        F = pic::getMatrixdFromLinearArray(F_data_opt, 3, 3);

        printf("Ok.\n");

        printf("\nFoundamental matrix: \n");
        pic::MatrixConvert(F).print();

        //Essential matrix decomposition
        Eigen::Matrix3d E = pic::computeEssentialMatrix(F, K);

        std::vector< Eigen::Vector2f > m0f, m1f;
        pic::filterInliers(m0, inliers, m0f);
        pic::filterInliers(m1, inliers, m1f);

        //Estimating R and t
        Eigen::Matrix3d R;
        Eigen::Vector3d t;
        pic::decomposeEssentialMatrixWithConfiguration(E, K, K, m0f, m1f, R, t);

        //Triangulation
        pic::Image imgOut0(1, img0.width, img0.height, 3);
        imgOut0.SetZero();

        pic::Image imgOut1(1, img1.width, img1.height, 3);
        imgOut1.SetZero();

        std::vector<Eigen::Vector3d> points_3d;

        Eigen::Matrix34d M0 = pic::getCameraMatrixIdentity(K);
        Eigen::Matrix34d M1 = pic::getCameraMatrix(K, R, t);

        pic::NelderMeadOptTriangulation nmTri(M0, M1);
        for(unsigned int i = 0; i < m0f.size(); i++) {
            //normalized coordinates
            Eigen::Vector3d p0 = Eigen::Vector3d(m0f[i][0], m0f[i][1], 1.0);
            Eigen::Vector3d p1 = Eigen::Vector3d(m1f[i][0], m1f[i][1], 1.0);

            //triangulation
            Eigen::Vector4d point = pic::triangulationHartleySturm(p0, p1, M0, M1);

            //refinement
            nmTri.update(m0f[i], m1f[i]);
            float tmpp[] = {float(point[0]), float(point[1]), float(point[2])};
            float out[3];
            nmTri.run(tmpp, 3, 1e-9f, 10000, &out[0]);

            points_3d.push_back(Eigen::Vector3d(out[0], out[1], out[2]));
        }

        pic::NelderMeadOptRadialDistortion nmRD(M0, M1, &m0f, &m1f, &points_3d);

        float lambda = 0.0f;
        float lambda_out;
        nmRD.run(&lambda, 1, 1e-12f, 100000, &lambda_out);
        printf("Radial distortion lambda: %f\n", lambda_out);

        double cx = 2592.0 / 2.0;
        double cy = 1728.0 / 2.0;
        for(unsigned int i = 0; i < m0f.size(); i++) {
            //first image
            Eigen::Vector2i proj0 = Projection(M0, points_3d[i], cx, cy, fx, fy, lambda_out);
            //Eigen::Vector2i proj0 = pic::cameraMatrixProject(M0, points_3d[i]);
            float *tmp;

            tmp = imgOut0(int(m0f[i][0]), int(m0f[i][1]));
            tmp[1] = 1.0f;

            tmp = imgOut0(proj0[0], proj0[1]);
            tmp[0] = 1.0f;

            //second image
            Eigen::Vector2i proj1 = Projection(M1, points_3d[i], cx, cy, fx, fy, lambda_out);
            //Eigen::Vector2i proj1 = pic::cameraMatrixProject(M1, points_3d[i]);

            tmp = imgOut1(int(m1f[i][0]), int(m1f[i][1]));
            tmp[1] = 1.0f;

            tmp = imgOut1(proj1[0], proj1[1]);
            tmp[0] = 1.0f;
        }

        imgOut0.Write("../data/output/simple_triangulation_reprojection_left.png");
        imgOut1.Write("../data/output/simple_triangulation_reprojection_right.png");

        //writing a PLY file
        FILE *file = fopen("../data/output/simple_triangulation_mesh.ply","w");

        if (file == NULL) {
            return 0;
        }

        fprintf(file,"ply\n");
        fprintf(file,"format ascii 1.0\n");
        fprintf(file,"element vertex %d\n", int(m0f.size()));

        fprintf(file,"property float x\n");
        fprintf(file,"property float y\n");
        fprintf(file,"property float z\n");

        fprintf(file,"property uchar red\n");
        fprintf(file,"property uchar green\n");
        fprintf(file,"property uchar blue\n");
        fprintf(file,"property uchar alpha\n");
        fprintf(file,"end_header\n");

        for(unsigned int i = 0; i < m0f.size(); i++) {
            fprintf(file, "%3.4f %3.4f %3.4f ", points_3d[i][0], points_3d[i][1], points_3d[i][2]);

            //writing color information
            unsigned char r, g, b;
            float *color = img0(int(m0f[i][0]), int(m0f[i][1]));
            r = int(color[0] * 255.0f);
            g = int(color[1] * 255.0f);
            b = int(color[2] * 255.0f);
            fprintf(file, " %d %d %d 255\n", r, g, b);
        }

        fclose(file);
    } else {
        printf("No there is at least an invalid file!\n");
    }

    return 0;
}
