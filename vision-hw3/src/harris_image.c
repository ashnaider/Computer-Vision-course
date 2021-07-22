#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}


float gaussian_1d_value(int x, float sigma)
{
    /* float gaussian = (1/(TWOPI*pow(sigma, 2))) * (exp(-(pow(x, 2))/(2*pow(sigma, 2)))); */
    /* return gaussian; */
    return (1.0 / (TWOPI * pow(sigma, 2))) * (exp(-(pow(x, 2)) / (2 * pow(sigma, 2))));
}


// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
    float sigma6 = 6 * sigma;
    int size = (int)sigma6 + 1;

    if (size % 2 == 0) {
        ++size;
    }

    image filter = make_image(size,1,1);

    float curr;
    /* for (int j = size/2; j >= -size/2; --j) { */
    /*     curr = gaussian_1d_value(j, sigma); */
    /*     set_pixel(filter, j, 0, 0, curr); */
    /* } */

    for(int j=0; j<size; j++){
        curr = gaussian_1d_value((size/2)-j, sigma);
        set_pixel(filter, j, 0, 0, curr);
    }

    l1_normalize(filter);

    return filter;
}

image transpose_image(image im)
{
    image res = make_image(im.h, im.w, im.c);
    float curr;

    for (int c = 0; c < im.c; ++c){
        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {
                curr = get_pixel(im, j, i, c);
                set_pixel(res, i, j, c, curr);
            }
        }
    }

    return res;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(0){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.

        image row_gauss_filter = make_1d_gaussian(sigma);
        image col_gauss_filter = transpose_image(row_gauss_filter);

        image im_row = convolve_image(im, row_gauss_filter, 1);
        image res = convolve_image(im_row, col_gauss_filter, 1);

        free_image(row_gauss_filter);
        free_image(col_gauss_filter);
        free_image(im_row);

        return res;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
    image Ix_filter = make_gx_filter();
    image Iy_filter = make_gy_filter();

    image Ix = convolve_image(im, Ix_filter, 0);
    image Iy = convolve_image(im, Iy_filter, 0);

    int ch = 0;
    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {

            float ch0_pxl = get_pixel(Ix, j, i, 0);
            float ch1_pxl = get_pixel(Iy, j, i, 0);

            float ch2_pxl = ch0_pxl * ch1_pxl;

            ch0_pxl = ch0_pxl * ch0_pxl;
            ch1_pxl = ch1_pxl * ch1_pxl;

            set_pixel(S, j, i, 0, ch0_pxl);
            set_pixel(S, j, i, 1, ch1_pxl);
            set_pixel(S, j, i, 2, ch2_pxl);
        }
    }

    image res = smooth_image(S, sigma);

    free_image(Ix_filter);
    free_image(Iy_filter);
    free_image(Ix);
    free_image(Iy);

    free_image(S);

    return res;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.

    float curr, a = 1, b, c, D, det, trace, x11, x22, x12, eigen1, eigen2;
    float alpha = 0.06;
    for (int i = 0; i < S.h; ++i) {
        for (int j = 0; j < S.w; ++j) {
            x11 = get_pixel(S, j, i, 0);
            x22 = get_pixel(S, j, i, 1);
            x12 = get_pixel(S, j, i, 2);

            b = -(x11 + x22);
            c = x11 * x22 - x12 * x12;

            D = sqrt(b*b - 4 * a * c);
            eigen1 = (-b + D) / (2 * a);
            eigen2 = (-b - D) / (2 * a);

            det = eigen1 * eigen2;
            trace = eigen1 + eigen2;
            curr = det - alpha * trace * trace;

            set_pixel(R, j, i, 0, curr);

        }
    }

    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])

    float center, curr;
    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            center = get_pixel(im, j, i, 0);

            for (int k = i - w; k < i + w; ++k) {
                for (int l = j - w; l < j + w; ++l) {

                    if (k == i && l == j) {
                        continue;
                    } else {
                        curr = get_pixel(im, l, k, 0);
                        if (curr > center) {
                            set_pixel(r, j, i, 0, -999999);
                            goto done_with_curr_pixel;
                        }
                    }

                }
            }
            done_with_curr_pixel:;
        }
    }

    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0; // change this
    float curr;
    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            curr = get_pixel(Rnms, j, i, 0);
            if (curr > thresh) {
                ++count;
            }
        }
    }


    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.

    int index = 0;
    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            curr = get_pixel(Rnms, j, i, 0);
            if (curr > thresh) {
                d[index] = describe_index(im, i * im.w + j);
                ++index;
            }
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
