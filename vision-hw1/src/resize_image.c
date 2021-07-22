#include "test.h"
#include <math.h>
#include "image.h"

typedef struct Pixel {
    float x;
    float y;
} coords;

typedef struct Coefs {
    float a;
    float b;
} coeffs;

float nn_interpolate(image im, float x, float y, int c)
{
    int new_x = (int)x, new_y = (int)y;

    if (x - (int)x >= 0.5) {
        new_x += 1;
    }

    if (y - (int)y >= 0.5) {
        new_y += 1;
    }

    return get_pixel(im, new_x, new_y, c);
}


coeffs get_coefs(float new_coord_start, float new_coord_end, float old_coord_start, float old_coord_end) {
    coeffs res;
    res.a = (old_coord_end - old_coord_start) / (new_coord_end - new_coord_start);
    res.b = old_coord_start - res.a * new_coord_start;
    return res;
}

coords map2old(coords curr, coeffs h_coef, coeffs w_coef) {
    coords res;
    res.x = w_coef.a * curr.x + w_coef.b;
    res.y = h_coef.a * curr.y + h_coef.b;
    return res;
}


image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image res = make_image(w, h, im.c);
    float new_val;
    coords old_pos, curr_pos;

    coeffs h_coef = get_coefs(-0.5, h-0.5, -0.5, im.h-0.5);
    coeffs w_coef = get_coefs(-0.5, w-0.5, -0.5, im.w-0.5);

    for (int ch = 0; ch < im.c; ++ch) {
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {

                curr_pos.x = j;
                curr_pos.y = i;

                old_pos = map2old(curr_pos, h_coef, w_coef);

                new_val = nn_interpolate(im, old_pos.x, old_pos.y, ch);

                set_pixel(res, j, i, ch, new_val);
            }
        }
    }


    return res;
}

float bilinear_interpolate(image im, float x, float y, int c)
{

    float d1, d2;

    d1 = y - (int)y;
    d2 = 1 - d1;

    float a_val, b_val, d_val, e_val;

    a_val = get_pixel(im, x, y, c);
    b_val = get_pixel(im, x + 1, y, c);
    d_val = get_pixel(im, x + 1, y + 1, c);
    e_val = get_pixel(im, x, y + 1, c);

    float q, q1, q2;

    q1 = d1 * e_val + d2 * a_val;
    q2 = d1 * d_val + d2 * b_val;


    d1 = x - (int)x;
    d2 = 1 - d1;

    q = d1 * q2 + d2 * q1;

    return q;
}

image bilinear_resize(image im, int w, int h)
{
    image res = make_image(w, h, im.c);
    float new_val;
    coords old_pos, curr_pos;

    coeffs h_coef = get_coefs(-0.5, h-0.5, -0.5, im.h-0.5);
    coeffs w_coef = get_coefs(-0.5, w-0.5, -0.5, im.w-0.5);

    for (int ch = 0; ch < im.c; ++ch) {
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {

                curr_pos.x = j;
                curr_pos.y = i;

                old_pos = map2old(curr_pos, h_coef, w_coef);

                new_val = bilinear_interpolate(im, old_pos.x, old_pos.y, ch);

                set_pixel(res, j, i, ch, new_val);

            }
        }
    }

    return res;
}

