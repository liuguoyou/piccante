/*

PICCANTE
The hottest HDR imaging library!
http://vcg.isti.cnr.it/piccante

Copyright (C) 2014
Visual Computing Laboratory - ISTI CNR
http://vcg.isti.cnr.it
First author: Francesco Banterle

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

*/

#ifndef PIC_GL_FILTERING_REINHARD_TMO_SINGLE_PASS_HPP
#define PIC_GL_FILTERING_REINHARD_TMO_SINGLE_PASS_HPP

#include "gl/filtering/filter.hpp"
#include "util/file_lister.hpp"
#include "gl/point_samplers/sampler_random_m.hpp"

namespace pic {

/**
 * @brief The FilterGLReinhardSinglePass class
 */
class FilterGLReinhardSinglePass: public FilterGL
{
protected:
    float sigma_s, sigma_r;
    MRSamplersGL<2> *ms;

    //tmo
    float alpha;

    //Random numbers tile
    ImageGL *imageRand;

    void InitShaders();
    void FragmentShader();

public:
    float Lwa;
    /**
     * @brief FilterGLReinhardSinglePass
     * @param sigma_s
     * @param sigma_r
     * @param type
     */
    FilterGLReinhardSinglePass(float alpha, float phi);

    ~FilterGLReinhardSinglePass();

    /**
     * @brief Update
     * @param sigma_s
     * @param sigma_r
     */
    void Update(float sigma_s, float sigma_r, float Lwa);

    /**
     * @brief Process
     * @param imgIn
     * @param imgOut
     * @return
     */
    ImageGL *Process(ImageGLVec imgIn, ImageGL *imgOut);
\

};

FilterGLReinhardSinglePass::FilterGLReinhardSinglePass(float alpha, float phi = 8.0f): FilterGL()
{
    this->alpha = alpha;

    float epsilon = 0.05f;
    float s_max = 8.0f;
    float sigma_s = 0.56f * powf(1.6f, s_max);
    float sigma_r = (powf(2.0f, phi) * alpha / (s_max * s_max)) * epsilon;

    //protected values are assigned/computed
    this->sigma_s = sigma_s;
    this->sigma_r = sigma_r;

    //Precomputation of the Gaussian Kernel
    int kernelSize = PrecomputedGaussian::KernelSize(sigma_s);//,sigma_r);
    int halfKernelSize = kernelSize >> 1;

    //Random numbers
    int nSamplers;

    imageRand = new ImageGL(1, 128, 128, 1, IMG_CPU, GL_TEXTURE_2D);
    imageRand->SetRand();
    imageRand->loadFromMemory();
    *imageRand -= 0.5f;
    nSamplers = 1;


    //Poisson samples
#ifdef PIC_DEBUG
    printf("Window: %d\n", halfKernelSize);
#endif

    ms = new MRSamplersGL<2>(ST_BRIDSON, halfKernelSize, halfKernelSize, 1,
                             nSamplers);
    ms->generateTexture();

    FragmentShader();
    InitShaders();
}

FilterGLReinhardSinglePass::~FilterGLReinhardSinglePass()
{
    delete imageRand;
    delete ms;

    //free shader etc...
}

void FilterGLReinhardSinglePass::FragmentShader()
{
    fragment_source = GLW_STRINGFY
                                          (
                                                  uniform sampler2D  u_tex;
                                                  uniform sampler2D  u_tex_col;
                                                  uniform isampler2D u_poisson;
                                                  uniform sampler2D  u_rand;
                                                  uniform int   nSamples;
                                                  uniform float sigmas2;
                                                  uniform float sigmar2;
                                                  uniform int kernelSize;
                                                  uniform float kernelSizef;
                                                  uniform float alpha;
                                                  out     vec4  f_color;

    void main(void) {
        ivec2 coordsFrag = ivec2(gl_FragCoord.xy);

        float colRef = texelFetch(u_tex, coordsFrag, 0).x;
        float Lw = colRef;

        colRef = colRef / (colRef + 1.0);


        float shifter = texture(u_rand, gl_FragCoord.xy).x;
        float  color = 0.0;
        float weight = 0.0;

        for(int i = 0; i < nSamples; i++) {
            //Coordinates
            ivec3 coords = texelFetch(u_poisson, ivec2(i, shifter), 0).xyz;


            //Texture fetch
            float tmpCol = texelFetch(u_tex, coordsFrag.xy + coords.xy, 0).x;
            tmpCol = tmpCol / (tmpCol + 1.0);

            float tmpCol2 = tmpCol - colRef;
            float dstR = tmpCol2 * tmpCol2;

            int coordsz = coords.x * coords.x + coords.y * coords.y;
            float tmp = exp(-dstR / sigmar2 - float(coordsz) / sigmas2);

            color += tmpCol * tmp;
            weight += tmp;

        }


        float bilateral = weight > 0.0 ? (color / weight) : colRef;
        bilateral = bilateral / (1.0 - bilateral);


        vec3 color_hdr = texelFetch(u_tex_col, coordsFrag, 0).xyz;

        float Ld = (Lw * alpha) / (bilateral * alpha + 1.0);

        f_color = vec4(color_hdr * (Ld / Lw), 1.0);


    }
                                          );

}

void FilterGLReinhardSinglePass::InitShaders()
{
#ifdef PIC_DEBUG
    printf("Number of samples: %d\n", ms->nSamples);
#endif

    filteringProgram.setup(glw::version("330"), vertex_source, fragment_source);

#ifdef PIC_DEBUG
    printf("[filteringProgram log]\n%s\n", filteringProgram.log().c_str());
#endif

    glw::bind_program(filteringProgram);
    filteringProgram.attribute_source("a_position", 0);
    filteringProgram.fragment_target("f_color",    0);
    filteringProgram.relink();
    glw::bind_program(0);

    Update(-1.0f, -1.0f, 1.0f);
}

void FilterGLReinhardSinglePass::Update(float sigma_s, float sigma_r, float Lwa)
{
    bool flag = false;
    this->Lwa = Lwa;

    if(sigma_s > 0.0f) {
        flag = (this->sigma_s == sigma_s);
        this->sigma_s = sigma_s;
    }

    if(sigma_r > 0.0f) {
        flag = flag || (this->sigma_r == sigma_r);
        this->sigma_r = sigma_r;
    }

    int kernelSize = PrecomputedGaussian::KernelSize(this->sigma_s);
    int halfKernelSize = kernelSize >> 1;

    if(flag) {
        ms->updateGL(halfKernelSize, halfKernelSize);
    }

    //shader update
    float sigmas2 = 2.0f * this->sigma_s * this->sigma_s;
    float sigmar2 = 2.0f * this->sigma_r * this->sigma_r;

    glw::bind_program(filteringProgram);
    filteringProgram.uniform("u_tex",       0);
    filteringProgram.uniform("u_poisson",   1);
    filteringProgram.uniform("u_rand",      2);
    filteringProgram.uniform("u_tex_col",      3);

    filteringProgram.uniform("sigmas2",         sigmas2);
    filteringProgram.uniform("alpha",           alpha / Lwa);
    filteringProgram.uniform("sigmar2",         sigmar2);
    filteringProgram.uniform("kernelSize",      kernelSize);
    filteringProgram.uniform("kernelSizef",     float(kernelSize));
    filteringProgram.uniform("nSamples",        ms->nSamples >> 1);
    glw::bind_program(0);
}

//Processing
ImageGL *FilterGLReinhardSinglePass::Process(ImageGLVec imgIn,
        ImageGL *imgOut)
{
    if(imgIn[0] == NULL) {
        return imgOut;
    }

    int w = imgIn[0]->width;
    int h = imgIn[0]->height;

    //TODO: check if other have height and frames swapped
    if(imgOut == NULL) {
        imgOut = new ImageGL(imgIn[0]->frames, w, h, imgIn[0]->channels, IMG_GPU, GL_TEXTURE_2D);
    }

    if(fbo == NULL) {
        fbo = new Fbo();
    }

    fbo->create(w, h, imgIn[0]->frames, false, imgOut->getTexture());

    ImageGL *edge, *base;

    edge = imgIn[0];
    base = imgIn[1];

    //Rendering
    fbo->bind();
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);

    //Shaders
    glw::bind_program(filteringProgram);

    //Textures
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, edge->getTexture());

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, imageRand->getTexture());

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, ms->getTexture());

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, base->getTexture());

    //Rendering aligned quad
    quad->Render();

    //Fbo
    fbo->unbind();

    //Shaders
    glw::bind_program(0);

    //Textures
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, 0);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, 0);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, 0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, 0);

    return imgOut;
}

} // end namespace pic

#endif /* PIC_GL_FILTERING_REINHARD_TMO_SINGLE_PASS_HPP */

