#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

int calculate_coord(image im, int x, int y, int c) {
    return c * im.h * im.w + y * im.w + x;
}

float get_pixel(image im, int x, int y, int c)
{
    if (x < 0) {
        x = 0;
    }
    if (x >= im.w) {
        x = im.w - 1;
    }
    if (y < 0) {
        y = 0;
    }
    if (y >= im.h) {
        y = im.h - 1;
    }

    int coord = calculate_coord(im,x, y, c);
    return im.data[coord];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (x < 0 || x >= im.w ||
        y < 0 || y >= im.h) {
        return ;
    }

    int coord = calculate_coord(im, x, y, c);
    im.data[coord] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, sizeof(float) * im.w * im.h * im.c);
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3); // what it does?
    image gray = make_image(im.w, im.h, 1);

    float red_comp, green_comp, blue_comp, res;

    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {

            red_comp = get_pixel(im, j, i, 0) * 0.299;
            green_comp = get_pixel(im, j, i, 1) * 0.587;
            blue_comp = get_pixel(im, j, i, 2) * 0.114;

            res = red_comp + green_comp + blue_comp;

            set_pixel(gray, j, i, 0, res);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    float curr;

    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {

            curr = get_pixel(im, j, i, c);
            set_pixel(im, j, i, c, curr + v);

        }
    }

}

void scale_image(image im, int c, float v) {
    float curr;

    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            curr = get_pixel(im, j, i, c);
            curr *= v;
            set_pixel(im, j, i, c, curr);
        }
    }
}

void clamp_image(image im)
{

    for (int c = 0; c < im.c; ++c) {
        for (int i = 0; i < im.h; ++i) {
            for (int j = 0; j < im.w; ++j) {

                if (get_pixel(im, j, i, c) > 1) {
                    set_pixel(im, j, i, c, 1);
                }
                if (get_pixel(im, j, i, c) < 0) {
                    set_pixel(im, j, i, c, 0);
                }
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    float V, m, C, S, dH, H;

    float R, G, B;

    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {

            R = get_pixel(im, j, i, 0);
            G = get_pixel(im, j, i, 1);
            B = get_pixel(im, j, i, 2);

            V = three_way_max(R, G, B);

            if (R == 0 && G == 0 && B == 0) {
                S = 0;
                H = 0;

                set_pixel(im, j, i, 0, H);
                set_pixel(im, j, i, 1, S);
                set_pixel(im, j, i, 2, V);

                continue;

            } else {
                m = three_way_min(R, G, B);

                C = V - m;
                S = C / V;
            }

            if (C == 0 || (V == 0 && m == 0)) {
                H = 0;

                set_pixel(im, j, i, 0, H);
                set_pixel(im, j, i, 1, S);
                set_pixel(im, j, i, 2, V);

                continue;
            }

            if (V == R) {
                dH = (G - B) / C;
            }
            else if (V == G) {
                dH = (B - R) / C + 2;
            }
            else if (V == B) {
                dH = (R - G) / C + 4;
            }



            if (dH < 0) {
                H = dH / 6 + 1;
            } else {
                H = dH / 6;
            }

            set_pixel(im, j, i, 0, H);
            set_pixel(im, j, i, 1, S);
            set_pixel(im, j, i, 2, V);

        }
    }
}

void hsv_to_rgb(image im)
{

    float C, m, X,
        r1, g1, b1,
        h, s, v,
        r, g, b;

    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {

            h = get_pixel(im, j, i, 0);
            s = get_pixel(im, j, i, 1);
            v = get_pixel(im, j, i, 2);

            C = v * s;

            h *= 6.0;

            X = C * (1 - fabs(fmod(h, 2) - 1));

            if (0 <= h && h <= 1.0) {
                r1 = C;
                g1 = X;
                b1 = 0;
            }
            else if (1.0 <= h && h <= 2.0) {
                r1 = X;
                g1 = C;
                b1 = 0;
            }
            else if (2.0 <= h && h <= 3.0) {
                r1 = 0;
                g1 = C;
                b1 = X;
            }
            else if (3.0 <= h && h <= 4.0) {
                r1 = 0;
                g1 = X;
                b1 = C;
            }
            else if (4.0 <= h && h <= 5.0) {
                r1 = X;
                g1 = 0;
                b1 = C;
            }
            else if (5.0 <= h && h <= 6.0) {
                r1 = C;
                g1 = 0;
                b1 = X;
            } else {
                r1 = 0;
                g1 = 0;
                b1 = 0;
            }

            m = v - C;

            r = r1 + m;
            g = g1 + m;
            b = b1 + m;

            set_pixel(im, j, i, 0, r);
            set_pixel(im, j, i, 1, g);
            set_pixel(im, j, i, 2, b);
        }
    }
}
