#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "test.h"
#define TWOPI 6.2831853

typedef struct MaxMin {
    float min;
    float max;
} min_max;

void l1_normalize(image im)
{
    float curr, sum = 0; // float square = im.w * im.h;

    for (int ch = 0; ch < im.c; ++ch) {
        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {
                sum += get_pixel(im, j, i, ch);
            }
        }
    }

    for (int c = 0; c < im.c; ++c) {
        scale_image(im, c, 1.0/sum);
    }
}

void l2_normalize(image im, float sum)
{
    float curr;
    for (int ch = 0; ch < im.c; ++ch) {
        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.h; ++j) {
                curr = get_pixel(im, j, i, ch);
                curr /= sum;
                set_pixel(im, j, i, ch, curr);
            }
        }
    }
}

image make_box_filter(int w)
{
    image filter = make_image(w, w, 1);

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < w; ++j) {
            set_pixel(filter, j, i, 0, 1);
        }
    }

    l1_normalize(filter);

    return filter;
}

float calc_filter(image im, int i, int j, int im_c, image filter, int filter_c) {
    float sum = 0;
    for (int fi = 0; fi < filter.h; ++fi) {
        for (int fj = 0; fj < filter.w; ++fj) {
            sum += get_pixel(filter, fj, fi, filter_c) * get_pixel(im, j - filter.w / 2 + fj, i - filter.h / 2 + fi, im_c);
        }
    }
    return sum;
}

image convolve_image(image im, image filter, int preserve)
{
    assert((im.c == filter.c) || filter.c == 1);

    image res;

    if (preserve == 1) {
        /* printf("in first case\n"); */
        res = make_image(im.w, im.h, im.c);
        float new_val, filter_ch;

        for (int ch = 0; ch < im.c; ++ch) {
            for (int i = 0; i < im.h; ++i) {
                for (int j = 0; j < im.w; ++j) {

                    if (filter.c == 1) {
                        filter_ch = 0;
                    } else {
                        filter_ch = ch;
                    }

                    new_val = calc_filter(im, i, j, ch, filter, filter_ch);
                    set_pixel(res, j, i, ch, new_val);

                }
            }
        }
    }
    else if (filter.c == 1 && im.c > 1) {
        res = make_image(im.w, im.h, 1);

        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {

                float sum = 0;

                for (int ch = 0; ch < im.c; ++ch) {
                    sum += calc_filter(im, i, j, ch, filter, 0);
                }

                set_pixel(res, j, i, 0, sum);

            } 
        }
    }

    else if (im.c == filter.c) {
        res = make_image(im.w, im.h, 1);
        float sum;

        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {

                sum = 0;

                for (int ch = 0; ch < im.c; ++ch) {
                    sum += calc_filter(im, i, j, ch, filter, ch);

                }

                set_pixel(res, j, i, 0, sum);
            }
        }
    }

    return res;
}

image make_highpass_filter()
{
    image filter = make_image(3, 3, 1);
    int ch = 0;
    set_pixel(filter, 0, 0, ch, 0);
    set_pixel(filter, 1, 0, ch, -1);
    set_pixel(filter, 2, 0, ch, 0);

    set_pixel(filter, 0, 1, ch, -1);
    set_pixel(filter, 1, 1, ch, 4);
    set_pixel(filter, 2, 1, ch, -1);

    set_pixel(filter, 0, 2, ch, 0);
    set_pixel(filter, 1, 2, ch, -1);
    set_pixel(filter, 2, 2, ch, 0);

    return filter;
}

image make_sharpen_filter()
{
    image filter = make_image(3, 3, 1);
    int ch = 0;
    set_pixel(filter, 0, 0, ch, 0);
    set_pixel(filter, 1, 0, ch, -1);
    set_pixel(filter, 2, 0, ch, 0);

    set_pixel(filter, 0, 1, ch, -1);
    set_pixel(filter, 1, 1, ch, 5);
    set_pixel(filter, 2, 1, ch, -1);

    set_pixel(filter, 0, 2, ch, 0);
    set_pixel(filter, 1, 2, ch, -1);
    set_pixel(filter, 2, 2, ch, 0);

    return filter;
}

image make_emboss_filter()
{
    image filter = make_image(3, 3, 1);
    int ch = 0;
    set_pixel(filter, 0, 0, ch, -2);
    set_pixel(filter, 1, 0, ch, -1);
    set_pixel(filter, 2, 0, ch, 0);

    set_pixel(filter, 0, 1, ch, -1);
    set_pixel(filter, 1, 1, ch, 1);
    set_pixel(filter, 2, 1, ch, 1);

    set_pixel(filter, 0, 2, ch, 0);
    set_pixel(filter, 1, 2, ch, 1);
    set_pixel(filter, 2, 2, ch, 2);

    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

float gaussian_value(int x, int y, float sigma) {
    float sigma_square = pow(sigma, 2);
    return exp(-1.0 * ((pow(x, 2) + pow(y, 2)) / (2.0 * sigma_square))) / (TWOPI * sigma_square);
}

image make_gaussian_filter(float sigma)
{
    float sigma6 = 6 * sigma;
    int size = (int)sigma6 + 1;

    if (size % 2 == 0) {
        ++size;
    }

    image filter = make_image(size, size, 1);
    float sum = 0, curr;

    int b = (size / 2);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            curr = gaussian_value(j-b, i-b, sigma);
            set_pixel(filter, j, i, 0, curr);

            sum += curr;
        }
    }

    l2_normalize(filter, sum);

    return filter;
}

image add_image(image a, image b)
{
    assert((a.w == b.w) && (a.h == b.h));

    assert(a.c == b.c);

    image res = make_image(a.w, a.h, a.c);
    float curr;

    for (int ch = 0; ch < a.c; ++ch) {
        for (int i = 0; i < a.h; ++i) {
            for (int j = 0; j < a.w; ++j) {
                curr = get_pixel(a, j, i, ch) + get_pixel(b, j, i, ch);

                set_pixel(res, j, i, ch, curr);
            }
        }
    }

    return res;
}

image sub_image(image a, image b)
{
    assert((a.w == b.w) && (a.h == b.h));

    assert(a.c == b.c);

    image res = make_image(a.w, a.h, a.c);
    float curr;

    for (int ch = 0; ch < a.c; ++ch) {
        for (int i = 0; i < a.h; ++i) {
            for (int j = 0; j < a.w; ++j) {
                curr = get_pixel(a, j, i, ch) - get_pixel(b, j, i, ch);

                set_pixel(res, j, i, ch, curr);
            }
        }
    }

    return res;
}

image make_gx_filter()
{
    image filter = make_image(3, 3, 1);
    int ch = 0;
    set_pixel(filter, 0, 0, ch, -1);
    set_pixel(filter, 1, 0, ch, 0);
    set_pixel(filter, 2, 0, ch, 1);

    set_pixel(filter, 0, 1, ch, -2);
    set_pixel(filter, 1, 1, ch, 0);
    set_pixel(filter, 2, 1, ch, 2);

    set_pixel(filter, 0, 2, ch, -1);
    set_pixel(filter, 1, 2, ch, 0);
    set_pixel(filter, 2, 2, ch, 1);

    return filter;
}

image make_gy_filter()
{
    image filter = make_image(3, 3, 1);
    int ch = 0;
    set_pixel(filter, 0, 0, ch, -1);
    set_pixel(filter, 1, 0, ch, -2);
    set_pixel(filter, 2, 0, ch, -1);

    set_pixel(filter, 0, 1, ch, 0);
    set_pixel(filter, 1, 1, ch, 0);
    set_pixel(filter, 2, 1, ch, 0);

    set_pixel(filter, 0, 2, ch, 1);
    set_pixel(filter, 1, 2, ch, 2);
    set_pixel(filter, 2, 2, ch, 1);

    return filter;
}

min_max find_min_max(image im) {
    min_max res;
    res.min = 2;
    res.max = -1;
    float curr;

    for (int ch = 0; ch < im.c; ++ch) {
        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {
                curr = get_pixel(im, j, i, ch);
                if (curr < res.min)
                    res.min = curr;
                if (curr > res.max)
                    res.max = curr;
            }
        }
    }
    return res;
}


void feature_normalize(image im)
{
    min_max mm = find_min_max(im);
    float range = mm.max - mm.min, curr;
    if (within_eps(range, 0)) {
        range = 1;
    }

    for (int ch = 0; ch < im.c; ++ch) {
        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {
                curr = get_pixel(im, j, i, ch);
                curr -= mm.min;
                curr /= range;
                set_pixel(im, j, i, ch, curr);
            }
        }
    }
}

image *sobel_image(image im)
{
    image magnitude = make_image(im.w, im.h, 1),
        direction = make_image(im.w, im.h, 1);

    image gx_filter = make_gx_filter(),
        gy_filter = make_gy_filter();

    image gx_res = convolve_image(im, gx_filter, 0),
        gy_res = convolve_image(im, gy_filter, 0);

    float gx, gy, g, dir;
    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            gx = get_pixel(gx_res, j, i, 0);
            gy = get_pixel(gy_res, j, i, 0);

            g = sqrt(pow(gx, 2) + pow(gy, 2));

            dir = atan2(gy, gx);

            set_pixel(magnitude, j, i, 0, g);
            set_pixel(direction, j, i, 0, dir);
        }
    }

    free_image(gx_res);
    free_image(gy_res);

    free_image(gx_filter);
    free_image(gy_filter);

    image* res = calloc(2, sizeof(image));
    res[0] = magnitude;
    res[1] = direction;

    return res;
}

image colorize_sobel(image im)
{
    image res = copy_image(im);

    rgb_to_hsv(res);

    image* sobel = sobel_image(im);



    float mag, angle;
    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            mag = get_pixel(sobel[0], j, i, 0);
            angle = get_pixel(sobel[1], j, i, 0);

            set_pixel(res, j, i, 0, angle);
            set_pixel(res, j, i, 1, mag);
        } 
    }

    feature_normalize(sobel[0]);
    feature_normalize(sobel[1]);

    free(sobel);

    hsv_to_rgb(res);

    image gauss_filter = make_gaussian_filter(2);
    convolve_image(res, gauss_filter, 1);

    free_image(gauss_filter);

    return res;
}
