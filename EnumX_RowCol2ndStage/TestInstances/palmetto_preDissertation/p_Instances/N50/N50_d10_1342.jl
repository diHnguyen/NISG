global arcs = [1 9; 1 20; 1 30; 2 15; 2 17; 2 29; 2 38; 2 46; 3 4; 3 37; 3 47; 4 12; 4 18; 4 24; 4 27; 4 28; 4 29; 4 45; 5 9; 5 17; 5 21; 5 29; 5 31; 6 11; 6 14; 6 39; 7 4; 7 16; 7 17; 7 29; 8 3; 8 12; 8 26; 8 29; 8 46; 9 2; 9 19; 9 23; 9 27; 10 2; 10 6; 10 7; 10 9; 10 13; 10 21; 10 23; 10 27; 10 46; 10 47; 11 12; 11 14; 11 19; 11 30; 11 33; 11 35; 11 42; 12 22; 12 42; 12 48; 12 50; 13 2; 13 3; 13 11; 13 22; 13 23; 13 40; 13 41; 13 43; 13 48; 14 5; 14 23; 15 13; 15 20; 15 26; 15 33; 15 40; 16 15; 16 28; 16 41; 17 4; 17 7; 17 14; 17 20; 17 25; 17 38; 17 41; 17 43; 17 50; 18 8; 18 21; 19 3; 19 7; 19 24; 19 31; 19 49; 20 17; 20 25; 20 42; 20 43; 20 44; 21 24; 21 28; 21 46; 22 3; 22 17; 22 41; 23 6; 23 25; 23 32; 23 36; 23 45; 23 47; 24 10; 24 12; 24 18; 24 19; 24 33; 24 50; 25 7; 25 19; 25 43; 25 46; 26 19; 26 33; 26 36; 26 44; 27 16; 27 26; 28 14; 28 26; 28 34; 29 4; 29 14; 29 31; 29 46; 30 5; 30 14; 30 18; 30 23; 30 28; 30 33; 30 39; 31 5; 31 11; 31 20; 32 5; 32 16; 32 24; 32 33; 32 46; 33 2; 33 7; 33 9; 33 13; 33 20; 33 27; 33 48; 33 50; 34 4; 34 12; 34 18; 34 21; 34 24; 34 35; 34 47; 34 48; 35 12; 35 29; 36 13; 36 39; 36 45; 36 50; 37 19; 37 29; 37 41; 37 44; 38 28; 38 29; 38 49; 39 4; 39 22; 39 25; 39 31; 39 44; 39 48; 39 49; 40 6; 40 14; 40 17; 40 19; 41 7; 41 8; 41 15; 41 42; 42 21; 42 36; 42 37; 43 4; 43 10; 43 16; 43 21; 43 38; 44 3; 44 15; 44 37; 44 40; 44 49; 45 22; 45 28; 46 12; 46 15; 46 29; 46 34; 46 50; 47 3; 47 20; 47 22; 47 41; 47 48; 48 12; 48 39; 48 40; 48 43; 48 45; 48 49; 48 50; 49 2; 49 3; 49 18; 49 19; 49 23; 49 44; 49 46]
global d_x = [4.0, 1.0, 6.0, 7.0, 10.0, 2.0, 9.0, 10.0, 4.0, 9.0, 1.0, 5.0, 7.0, 1.0, 5.0, 7.0, 4.0, 10.0, 3.0, 6.0, 1.0, 8.0, 6.0, 5.0, 9.0, 10.0, 9.0, 5.0, 5.0, 2.0, 9.0, 7.0, 10.0, 8.0, 3.0, 1.0, 5.0, 10.0, 2.0, 3.0, 5.0, 8.0, 2.0, 6.0, 2.0, 5.0, 1.0, 2.0, 7.0, 6.0, 10.0, 7.0, 10.0, 9.0, 8.0, 9.0, 6.0, 6.0, 2.0, 10.0, 1.0, 7.0, 2.0, 5.0, 7.0, 10.0, 10.0, 1.0, 9.0, 9.0, 7.0, 7.0, 1.0, 2.0, 5.0, 3.0, 8.0, 5.0, 10.0, 3.0, 8.0, 9.0, 2.0, 2.0, 9.0, 8.0, 8.0, 2.0, 8.0, 10.0, 5.0, 4.0, 4.0, 4.0, 4.0, 9.0, 1.0, 5.0, 8.0, 4.0, 6.0, 9.0, 6.0, 5.0, 8.0, 4.0, 3.0, 6.0, 8.0, 6.0, 1.0, 10.0, 5.0, 10.0, 6.0, 5.0, 6.0, 4.0, 6.0, 10.0, 9.0, 1.0, 8.0, 2.0, 1.0, 5.0, 7.0, 3.0, 6.0, 6.0, 3.0, 1.0, 4.0, 7.0, 2.0, 10.0, 9.0, 10.0, 8.0, 9.0, 7.0, 2.0, 1.0, 5.0, 2.0, 5.0, 4.0, 7.0, 10.0, 6.0, 3.0, 3.0, 10.0, 9.0, 6.0, 9.0, 10.0, 5.0, 1.0, 9.0, 4.0, 2.0, 4.0, 6.0, 6.0, 3.0, 4.0, 10.0, 10.0, 4.0, 5.0, 7.0, 10.0, 3.0, 5.0, 7.0, 5.0, 1.0, 3.0, 8.0, 6.0, 5.0, 1.0, 6.0, 2.0, 10.0, 1.0, 3.0, 5.0, 4.0, 10.0, 10.0, 9.0, 3.0, 10.0, 3.0, 4.0, 1.0, 5.0, 7.0, 2.0, 10.0, 4.0, 2.0, 6.0, 6.0, 1.0, 5.0, 8.0, 5.0, 7.0, 2.0, 3.0, 8.0, 3.0, 1.0, 1.0, 7.0, 1.0, 9.0, 2.0, 2.0, 5.0, 1.0, 7.0, 1.0, 8.0, 9.0, 7.0, 3.0, 5.0, 1.0, 5.0]
global b_x = 5
global d_y = [7.0, 9.0, 5.0, 9.0, 4.0, 5.0, 9.0, 6.0, 1.0, 1.0, 10.0, 8.0, 9.0, 1.0, 4.0, 6.0, 4.0, 8.0, 6.0, 4.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0, 9.0, 7.0, 9.0, 9.0, 7.0, 1.0, 8.0, 5.0, 2.0, 2.0, 9.0, 4.0, 8.0, 2.0, 1.0, 4.0, 7.0, 3.0, 5.0, 6.0, 9.0, 2.0, 9.0, 7.0, 6.0, 7.0, 6.0, 9.0, 2.0, 1.0, 10.0, 5.0, 4.0, 4.0, 2.0, 4.0, 3.0, 7.0, 1.0, 7.0, 5.0, 2.0, 6.0, 2.0, 4.0, 9.0, 2.0, 3.0, 5.0, 8.0, 1.0, 4.0, 8.0, 5.0, 4.0, 8.0, 9.0, 4.0, 7.0, 6.0, 1.0, 9.0, 9.0, 1.0, 8.0, 9.0, 10.0, 4.0, 8.0, 9.0, 9.0, 1.0, 4.0, 9.0, 4.0, 7.0, 9.0, 6.0, 4.0, 8.0, 4.0, 8.0, 7.0, 5.0, 7.0, 6.0, 9.0, 7.0, 10.0, 10.0, 7.0, 10.0, 10.0, 3.0, 7.0, 8.0, 9.0, 1.0, 9.0, 1.0, 5.0, 9.0, 5.0, 3.0, 5.0, 9.0, 8.0, 1.0, 4.0, 3.0, 3.0, 1.0, 7.0, 5.0, 3.0, 1.0, 8.0, 2.0, 5.0, 6.0, 3.0, 8.0, 7.0, 2.0, 9.0, 9.0, 2.0, 9.0, 7.0, 9.0, 5.0, 8.0, 8.0, 6.0, 4.0, 3.0, 4.0, 8.0, 10.0, 5.0, 2.0, 10.0, 9.0, 1.0, 8.0, 6.0, 7.0, 8.0, 1.0, 3.0, 7.0, 10.0, 6.0, 2.0, 1.0, 6.0, 7.0, 3.0, 10.0, 10.0, 1.0, 4.0, 8.0, 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 9.0, 2.0, 9.0, 4.0, 7.0, 7.0, 3.0, 10.0, 10.0, 7.0, 4.0, 1.0, 9.0, 7.0, 6.0, 9.0, 10.0, 5.0, 8.0, 3.0, 7.0, 5.0, 7.0, 10.0, 5.0, 1.0, 10.0, 9.0, 1.0, 9.0, 5.0, 7.0, 9.0, 1.0, 8.0, 6.0, 9.0, 5.0]
global b_y = 10
global p = [0.366, 0.865, 0.864, 0.473, 0.292, 0.471, 0.901, 0.622, 0.822, 0.387, 0.596, 0.364, 0.871, 0.819, 0.19, 0.169, 0.342, 0.051, 0.052, 0.717, 0.497, 0.252, 0.191, 0.806, 0.295, 0.088, 0.603, 0.167, 0.312, 0.615, 0.069, 0.784, 0.279, 0.274, 0.617, 0.785, 0.183, 0.619, 0.728, 0.713, 0.675, 0.638, 0.369, 0.665, 0.067, 0.73, 0.508, 0.376, 0.241, 0.571, 0.365, 0.257, 0.345, 0.433, 0.568, 0.619, 0.868, 0.734, 0.756, 0.726, 0.148, 0.06, 0.855, 0.218, 0.733, 0.951, 0.613, 0.597, 0.339, 0.58, 0.079, 0.855, 0.98, 0.629, 0.32, 0.737, 0.939, 0.306, 0.253, 0.421, 0.829, 0.029, 0.136, 0.664, 0.622, 0.327, 0.281, 0.736, 0.151, 0.735, 0.334, 0.235, 0.846, 0.637, 0.595, 0.696, 0.315, 0.704, 0.092, 0.788, 0.558, 0.922, 0.508, 0.386, 0.492, 0.987, 0.421, 0.906, 0.951, 0.768, 0.837, 0.098, 0.329, 0.185, 0.687, 0.406, 0.4, 0.048, 0.161, 0.841, 0.079, 0.341, 0.99, 0.184, 0.103, 0.861, 0.929, 0.467, 0.231, 0.401, 0.235, 0.986, 0.347, 0.05, 0.849, 0.404, 0.012, 0.718, 0.164, 0.655, 0.099, 0.389, 0.824, 0.395, 0.297, 0.961, 0.845, 0.418, 0.528, 0.798, 0.376, 0.069, 0.258, 0.893, 0.521, 0.157, 0.568, 0.946, 0.009, 0.556, 0.564, 0.835, 0.326, 0.944, 0.276, 0.84, 0.055, 0.049, 0.829, 0.924, 0.258, 0.288, 0.902, 0.089, 0.689, 0.676, 0.718, 0.362, 0.935, 0.304, 0.865, 0.605, 0.822, 0.881, 0.506, 0.017, 0.009, 0.434, 0.74, 0.998, 0.856, 0.412, 0.938, 0.124, 0.779, 0.946, 0.692, 0.994, 0.235, 0.81, 0.46, 0.808, 0.567, 0.354, 0.578, 0.696, 0.673, 0.992, 0.289, 0.39, 0.967, 0.063, 0.956, 0.44, 0.362, 0.923, 0.931, 0.249, 0.183, 0.727, 0.24, 0.022, 0.892, 0.668, 0.013, 0.648, 0.405, 0.318, 0.636, 0.713, 0.808, 0.027, 0.805]
global q = [0.429, 0.934, 0.998, 0.906, 0.318, 0.952, 0.991, 0.75, 0.827, 0.52, 0.632, 0.878, 0.975, 0.865, 0.574, 0.193, 0.635, 0.43, 0.657, 0.93, 0.625, 0.736, 0.55, 0.827, 0.769, 0.996, 0.696, 0.684, 0.565, 0.93, 0.232, 0.857, 0.841, 0.291, 0.846, 0.886, 0.549, 0.876, 0.772, 0.953, 0.75, 0.996, 0.61, 0.9, 0.454, 0.888, 0.866, 0.508, 0.766, 0.796, 0.688, 0.57, 0.671, 0.857, 0.72, 0.821, 0.894, 0.775, 0.894, 0.745, 0.38, 0.6, 0.972, 0.229, 0.977, 0.98, 0.808, 0.99, 0.885, 0.768, 0.339, 0.997, 0.998, 0.722, 0.453, 0.941, 0.95, 0.878, 0.382, 0.579, 0.964, 0.228, 0.807, 0.807, 0.959, 0.885, 0.665, 0.746, 0.728, 0.947, 0.694, 0.973, 0.903, 0.668, 0.87, 0.995, 0.879, 0.854, 0.163, 0.788, 0.729, 0.933, 0.568, 0.896, 0.598, 0.998, 0.586, 0.958, 0.965, 0.953, 0.904, 0.878, 0.367, 0.631, 0.824, 0.503, 0.839, 0.259, 0.528, 0.944, 0.714, 0.831, 0.992, 0.733, 0.504, 0.925, 0.963, 0.673, 0.446, 0.654, 0.348, 0.988, 0.753, 0.392, 0.917, 0.775, 0.75, 0.86, 0.205, 0.955, 0.966, 0.884, 0.986, 0.663, 0.805, 0.98, 0.887, 0.696, 0.657, 0.84, 0.763, 0.244, 0.601, 0.969, 0.531, 0.843, 0.768, 0.95, 0.857, 0.685, 0.733, 0.967, 0.63, 0.979, 0.422, 0.998, 0.403, 0.715, 0.852, 0.986, 0.972, 0.47, 0.939, 0.301, 0.992, 0.939, 0.967, 0.646, 0.967, 0.694, 0.898, 0.981, 0.954, 0.953, 0.54, 0.258, 0.162, 0.89, 0.81, 0.999, 0.921, 0.423, 0.945, 0.282, 0.86, 0.967, 0.938, 0.994, 0.359, 0.875, 0.816, 0.911, 0.81, 0.897, 0.878, 0.889, 0.994, 0.997, 0.769, 0.463, 0.972, 0.285, 0.99, 0.857, 0.776, 0.956, 0.947, 0.509, 0.718, 0.878, 0.68, 0.436, 0.997, 0.881, 0.041, 0.667, 0.872, 0.452, 0.746, 0.852, 0.932, 0.894, 0.921]
global origin = 1
global destination = 50