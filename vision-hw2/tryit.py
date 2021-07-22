from uwimg import *
# im = load_image("data/dogsmall.jpg")
# a = nn_resize(im, im.w*4, im.h*4)
# save_image(a, "dog4x-nn")
# b = nn_resize(im, 713, 467)
# save_image(b, "dog_custom_size")

# im = load_image("data/dogsmall.jpg")
# a = bilinear_resize(im, im.w*4, im.h*4)
# save_image(a, "dog4x-bl")

# im = load_image("data/dog.jpg")
# b = bilinear_resize(im, 713, 467)
# save_image(b, "dog_bl_custom_size")

# a = nn_resize(im, im.w//7, im.h//7)
# save_image(a, "dog7th-bl")


# im = load_image("data/dog.jpg")
# f = make_box_filter(7)
# blur = convolve_image(im, f, 1)
# clamp_image(blur)
# save_image(blur, "dog-box7")

# im = load_image("data/dog.jpg")
# f = make_box_filter(7)
# blur = convolve_image(im, f, 1)
# thumb = nn_resize(blur, blur.w//7, blur.h//7)
# save_image(thumb, "dogthumb")


# im = load_image("data/dog.jpg")
# f = make_gaussian_filter(2)
# blur = convolve_image(im, f, 1)
# clamp_image(blur)

# save_image(blur, "dog-gauss2")


# im = load_image("data/dog.jpg")
# f = make_highpass_filter()

# edges = convolve_image(im, f, 0);

# save_image(edges, "temp-dog")

# im = load_image("data/dog.jpg")
# f = make_gaussian_filter(2)
# lfreq = convolve_image(im, f, 1)
# hfreq = im - lfreq
# reconstruct = lfreq + hfreq
# save_image(lfreq, "low-frequency")
# save_image(hfreq, "high-frequency")
# save_image(reconstruct, "reconstruct")


# im = load_image("data/dog.jpg")
# res = sobel_image(im)
# mag = res[1]
# feature_normalize(mag)
# save_image(mag, "theta")


# im = load_image("data/dog.jpg")
# res = colorize_sobel(im)
# save_image(res, "colorized_sobel")


# im = load_image("data/Rainier1.png")
# S = structure_matrix(im, 2)
# R = cornerness_response(S)
# save_image(R, "cornerness_test")
# nms = nms_image(R, 3)
# # feature_normalize(nms)
# save_image(nms, "nms_test")


# im = load_image("data/Rainier1.png")
# detect_and_draw_corners(im, 2, 50, 3)
# save_image(im, "corners")

# a = load_image("data/Rainier1.png")
# b = load_image("data/Rainier2.png")
# m = find_and_draw_matches(a, b, 2, 50, 3)
# save_image(m, "matches")

im1 = load_image("data/Rainier1.png")
im2 = load_image("data/Rainier2.png")
pan = panorama_image(im1, im2, thresh=50)
print("after panorama")
save_image(pan, "easy_panorama")
