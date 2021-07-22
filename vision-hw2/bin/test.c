#include "../src/image.h"
#include "../src/matrix.h"
#include "../src/stb_image.h"
#include "../src/test.c"
#include "../src/args.h"
#include "../src/stb_image_write.h"
#include "../src/test.h"
#include <stdio.h>

int main() {
    printf("ok\n");

    image im = load_image("../figs/matches.jpg");
    free_image(im);

    return 0;
}
