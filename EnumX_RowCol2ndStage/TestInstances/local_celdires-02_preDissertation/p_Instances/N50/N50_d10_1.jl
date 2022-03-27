global arcs = [1 10; 1 14; 1 28; 1 31; 1 33; 1 40; 1 46; 1 47; 2 6; 2 10; 2 21; 2 22; 2 23; 2 33; 2 34; 2 50; 3 15; 3 25; 3 29; 3 35; 3 44; 3 46; 4 13; 4 16; 4 18; 4 21; 4 23; 4 31; 4 36; 5 6; 5 13; 5 27; 5 34; 5 38; 5 39; 6 3; 6 12; 6 24; 6 25; 6 27; 6 28; 6 34; 7 5; 7 9; 7 13; 7 19; 7 33; 8 25; 8 29; 8 31; 8 33; 8 42; 9 24; 10 24; 10 32; 10 33; 10 38; 10 50; 11 44; 12 2; 12 4; 12 5; 12 21; 12 25; 12 29; 12 46; 13 7; 13 17; 13 24; 13 46; 13 47; 14 40; 14 44; 15 11; 15 42; 15 45; 16 2; 16 4; 16 7; 16 12; 16 36; 16 47; 17 3; 17 19; 17 26; 17 34; 17 45; 17 46; 17 48; 18 2; 18 22; 18 40; 18 42; 19 28; 19 30; 19 39; 20 2; 20 5; 20 12; 20 19; 20 23; 20 34; 20 39; 20 46; 21 5; 21 37; 22 9; 22 28; 22 31; 23 2; 23 17; 23 18; 23 21; 23 22; 23 24; 23 28; 23 40; 24 2; 24 3; 24 4; 24 15; 24 16; 24 18; 24 37; 24 43; 24 46; 25 9; 25 19; 25 26; 25 29; 25 40; 25 46; 25 47; 26 21; 26 35; 26 44; 26 47; 26 50; 27 5; 27 19; 27 29; 27 34; 27 39; 28 2; 28 8; 28 13; 28 14; 28 31; 28 45; 29 6; 29 8; 29 23; 29 30; 29 39; 30 6; 30 11; 30 20; 30 22; 30 35; 30 48; 30 49; 31 2; 31 3; 31 20; 31 45; 31 49; 32 13; 32 18; 32 30; 32 31; 32 38; 32 46; 33 14; 33 43; 33 49; 34 4; 34 11; 34 13; 34 21; 34 38; 34 43; 34 45; 35 10; 35 12; 35 21; 35 22; 35 32; 35 38; 35 39; 35 41; 35 50; 36 18; 36 19; 36 30; 37 3; 37 13; 37 20; 37 29; 37 30; 37 34; 37 48; 37 49; 38 21; 38 24; 38 25; 38 31; 38 37; 38 41; 38 50; 39 14; 39 18; 39 35; 39 41; 39 50; 40 27; 41 2; 41 8; 41 9; 41 13; 41 42; 42 14; 42 20; 42 29; 42 33; 42 38; 43 4; 43 5; 43 18; 43 24; 43 34; 43 36; 43 40; 43 46; 43 48; 43 50; 44 10; 44 13; 44 25; 44 27; 44 28; 45 3; 45 19; 45 23; 45 30; 45 41; 45 43; 45 44; 45 48; 46 11; 46 15; 46 23; 46 39; 46 50; 47 14; 47 27; 47 32; 47 38; 48 9; 48 10; 48 25; 48 28; 48 34; 49 5; 49 6; 49 25; 49 32; 49 37; 49 50]
global d_x = [4.0, 3.0, 1.0, 4.0, 2.0, 1.0, 2.0, 1.0, 10.0, 1.0, 10.0, 9.0, 10.0, 5.0, 6.0, 8.0, 3.0, 9.0, 5.0, 4.0, 6.0, 9.0, 10.0, 8.0, 7.0, 8.0, 5.0, 7.0, 7.0, 5.0, 10.0, 10.0, 5.0, 3.0, 6.0, 7.0, 4.0, 10.0, 9.0, 8.0, 4.0, 9.0, 10.0, 6.0, 5.0, 8.0, 10.0, 8.0, 9.0, 10.0, 1.0, 2.0, 7.0, 2.0, 3.0, 7.0, 5.0, 5.0, 10.0, 6.0, 1.0, 9.0, 7.0, 9.0, 2.0, 8.0, 2.0, 2.0, 4.0, 3.0, 3.0, 8.0, 2.0, 2.0, 5.0, 6.0, 8.0, 9.0, 6.0, 5.0, 9.0, 10.0, 10.0, 1.0, 5.0, 9.0, 2.0, 4.0, 7.0, 8.0, 4.0, 4.0, 5.0, 2.0, 8.0, 9.0, 4.0, 10.0, 6.0, 6.0, 1.0, 5.0, 5.0, 1.0, 8.0, 6.0, 10.0, 5.0, 6.0, 8.0, 6.0, 7.0, 10.0, 3.0, 6.0, 6.0, 6.0, 5.0, 3.0, 7.0, 4.0, 4.0, 3.0, 4.0, 5.0, 5.0, 5.0, 2.0, 10.0, 10.0, 7.0, 4.0, 10.0, 9.0, 8.0, 6.0, 2.0, 4.0, 5.0, 4.0, 3.0, 7.0, 8.0, 2.0, 1.0, 6.0, 7.0, 5.0, 5.0, 7.0, 6.0, 6.0, 8.0, 2.0, 7.0, 9.0, 4.0, 5.0, 1.0, 4.0, 7.0, 4.0, 9.0, 8.0, 2.0, 1.0, 10.0, 10.0, 7.0, 1.0, 8.0, 4.0, 3.0, 1.0, 5.0, 7.0, 8.0, 8.0, 2.0, 6.0, 7.0, 6.0, 8.0, 10.0, 7.0, 3.0, 4.0, 10.0, 10.0, 4.0, 6.0, 8.0, 1.0, 4.0, 3.0, 3.0, 10.0, 5.0, 4.0, 4.0, 8.0, 6.0, 2.0, 6.0, 1.0, 1.0, 3.0, 2.0, 10.0, 3.0, 10.0, 7.0, 1.0, 3.0, 6.0, 8.0, 7.0, 1.0, 6.0, 2.0, 8.0, 10.0, 10.0, 9.0, 3.0, 3.0, 4.0, 7.0, 3.0, 8.0, 10.0, 4.0, 9.0, 8.0, 9.0, 5.0, 6.0, 4.0, 3.0, 7.0, 1.0, 6.0, 7.0, 2.0, 9.0, 9.0, 6.0, 10.0, 2.0, 10.0, 3.0, 3.0, 9.0, 9.0, 4.0, 5.0, 6.0, 4.0, 10.0, 7.0, 2.0, 8.0, 7.0, 5.0, 4.0, 1.0, 9.0, 9.0]
global b_x = 5
global d_y = [4.0, 3.0, 8.0, 3.0, 6.0, 2.0, 9.0, 4.0, 7.0, 9.0, 2.0, 6.0, 6.0, 10.0, 5.0, 3.0, 6.0, 9.0, 2.0, 5.0, 5.0, 10.0, 2.0, 9.0, 3.0, 10.0, 4.0, 3.0, 9.0, 10.0, 5.0, 8.0, 4.0, 4.0, 7.0, 6.0, 4.0, 10.0, 8.0, 1.0, 4.0, 2.0, 6.0, 2.0, 10.0, 7.0, 5.0, 5.0, 6.0, 8.0, 2.0, 2.0, 4.0, 10.0, 4.0, 10.0, 1.0, 2.0, 2.0, 10.0, 8.0, 2.0, 1.0, 1.0, 10.0, 5.0, 9.0, 10.0, 10.0, 4.0, 7.0, 2.0, 4.0, 6.0, 4.0, 9.0, 5.0, 10.0, 10.0, 6.0, 4.0, 5.0, 5.0, 7.0, 4.0, 5.0, 5.0, 3.0, 1.0, 9.0, 2.0, 3.0, 10.0, 6.0, 10.0, 4.0, 3.0, 6.0, 6.0, 2.0, 2.0, 10.0, 5.0, 3.0, 3.0, 2.0, 7.0, 5.0, 5.0, 3.0, 5.0, 2.0, 9.0, 8.0, 10.0, 6.0, 4.0, 2.0, 4.0, 1.0, 4.0, 8.0, 8.0, 5.0, 6.0, 2.0, 9.0, 10.0, 7.0, 2.0, 8.0, 5.0, 1.0, 10.0, 6.0, 2.0, 2.0, 10.0, 2.0, 3.0, 8.0, 3.0, 7.0, 3.0, 9.0, 4.0, 3.0, 10.0, 2.0, 8.0, 1.0, 9.0, 1.0, 5.0, 3.0, 1.0, 4.0, 6.0, 2.0, 8.0, 6.0, 9.0, 10.0, 4.0, 1.0, 9.0, 9.0, 2.0, 10.0, 10.0, 8.0, 5.0, 6.0, 6.0, 9.0, 7.0, 3.0, 7.0, 7.0, 4.0, 1.0, 1.0, 2.0, 7.0, 4.0, 4.0, 8.0, 7.0, 7.0, 8.0, 5.0, 6.0, 8.0, 7.0, 3.0, 7.0, 9.0, 4.0, 10.0, 3.0, 8.0, 7.0, 5.0, 8.0, 9.0, 5.0, 5.0, 5.0, 10.0, 2.0, 6.0, 4.0, 4.0, 5.0, 5.0, 5.0, 3.0, 8.0, 1.0, 1.0, 6.0, 2.0, 10.0, 2.0, 2.0, 9.0, 8.0, 1.0, 8.0, 10.0, 7.0, 9.0, 2.0, 10.0, 8.0, 6.0, 5.0, 10.0, 5.0, 8.0, 4.0, 7.0, 4.0, 5.0, 4.0, 4.0, 6.0, 9.0, 2.0, 7.0, 9.0, 10.0, 1.0, 6.0, 1.0, 7.0, 1.0, 10.0, 10.0, 2.0, 1.0, 6.0, 2.0, 10.0, 2.0, 1.0, 6.0, 6.0]
global b_y = 10
global p = [0.869, 0.708, 0.755, 0.961, 0.934, 0.843, 0.368, 0.993, 0.165, 0.225, 0.438, 0.806, 0.055, 0.39, 0.355, 0.909, 0.612, 0.535, 0.01, 0.246, 0.052, 0.79, 0.893, 0.263, 0.886, 0.536, 0.125, 0.712, 0.627, 0.037, 0.138, 0.996, 0.953, 0.828, 0.72, 0.541, 0.588, 0.989, 0.953, 0.839, 0.282, 0.135, 0.621, 0.463, 0.949, 0.981, 0.43, 0.259, 0.479, 0.239, 0.37, 0.461, 0.078, 0.19, 0.707, 0.597, 0.771, 0.928, 0.596, 0.847, 0.429, 0.855, 0.283, 0.424, 0.777, 0.165, 0.965, 0.064, 0.643, 0.853, 0.57, 0.012, 0.987, 0.182, 0.135, 0.831, 0.167, 0.576, 0.041, 0.606, 0.928, 0.722, 0.378, 0.236, 0.158, 0.03, 0.541, 0.375, 0.367, 0.465, 0.447, 0.528, 0.192, 0.186, 0.302, 0.405, 0.703, 0.347, 0.411, 0.293, 0.143, 0.186, 0.994, 0.186, 0.948, 0.257, 0.562, 0.305, 0.32, 0.892, 0.905, 0.437, 0.211, 0.151, 0.32, 0.444, 0.646, 0.671, 0.83, 0.078, 0.831, 0.153, 0.961, 0.085, 0.404, 0.156, 0.708, 0.154, 0.403, 0.278, 0.257, 0.276, 0.038, 0.276, 0.216, 0.794, 0.793, 0.186, 0.636, 0.308, 0.655, 0.972, 0.029, 0.253, 0.4, 0.575, 0.056, 0.832, 0.588, 0.839, 0.794, 0.254, 0.653, 0.945, 0.088, 0.069, 0.353, 0.278, 0.161, 0.368, 0.532, 0.009, 0.62, 0.22, 0.913, 0.471, 0.361, 0.712, 0.486, 0.713, 0.29, 0.924, 0.594, 0.804, 0.373, 0.381, 0.134, 0.819, 0.166, 0.241, 0.895, 0.556, 0.347, 0.785, 0.334, 0.653, 0.959, 0.017, 0.318, 0.925, 0.839, 0.443, 0.683, 0.386, 0.091, 0.254, 0.31, 0.992, 0.916, 0.273, 0.671, 0.27, 0.695, 0.334, 0.353, 0.553, 0.49, 0.522, 0.784, 0.103, 0.284, 0.815, 0.089, 0.116, 0.483, 0.619, 0.208, 0.844, 0.943, 0.24, 0.408, 0.84, 0.497, 0.873, 0.853, 0.001, 0.872, 0.289, 0.941, 0.391, 0.593, 0.919, 0.04, 0.83, 0.554, 0.032, 0.202, 0.528, 0.298, 0.317, 0.013, 0.565, 0.648, 0.967, 0.677, 0.85, 0.834, 0.314, 0.064, 0.007, 0.428, 0.679, 0.324, 0.176, 0.223, 0.008, 0.858, 0.403, 0.632, 0.892, 0.647, 0.293, 0.978, 0.953, 0.07, 0.221, 0.95, 0.841]
global q = [0.916, 0.745, 0.877, 0.97, 0.987, 0.981, 0.423, 0.999, 0.19, 0.586, 0.443, 0.91, 0.529, 0.828, 0.701, 0.926, 0.971, 0.668, 0.122, 0.482, 0.427, 0.901, 0.975, 0.405, 0.98, 0.97, 0.768, 0.852, 0.718, 0.886, 0.653, 0.996, 0.985, 0.9, 0.814, 0.689, 0.655, 0.998, 0.979, 0.964, 0.334, 0.46, 0.648, 0.847, 0.989, 0.983, 0.6, 0.66, 0.625, 0.553, 0.547, 0.706, 0.396, 0.566, 0.844, 0.714, 0.945, 0.949, 0.846, 0.853, 0.597, 0.932, 0.406, 0.802, 0.953, 0.454, 0.979, 0.329, 0.829, 0.868, 0.811, 0.68, 0.993, 0.575, 0.598, 0.931, 0.993, 0.968, 0.296, 0.982, 0.966, 0.862, 0.884, 0.418, 0.587, 0.053, 0.686, 0.858, 0.518, 0.939, 0.96, 0.796, 0.402, 0.661, 0.991, 0.654, 0.751, 0.548, 0.748, 0.678, 0.25, 0.761, 0.997, 0.3, 0.957, 0.287, 0.968, 0.972, 0.643, 0.899, 0.927, 0.765, 0.795, 0.816, 0.682, 0.898, 0.871, 0.903, 0.853, 0.088, 0.968, 0.798, 0.975, 0.305, 0.924, 0.864, 0.786, 0.449, 0.773, 0.948, 0.899, 0.939, 0.224, 0.322, 0.941, 0.924, 0.985, 0.447, 0.758, 0.455, 0.692, 0.986, 0.466, 0.866, 0.842, 0.599, 0.808, 0.881, 0.707, 0.986, 0.962, 0.349, 0.858, 0.974, 0.333, 0.604, 0.365, 0.45, 0.393, 0.613, 0.685, 0.083, 0.7, 0.386, 0.958, 0.67, 0.403, 0.824, 0.708, 0.729, 0.48, 0.993, 0.964, 0.847, 0.686, 0.599, 0.55, 0.926, 0.927, 0.457, 0.988, 0.986, 0.951, 0.966, 0.989, 0.658, 0.963, 0.639, 0.478, 0.959, 0.928, 0.514, 0.825, 0.859, 0.294, 0.769, 0.959, 0.992, 0.935, 0.712, 0.785, 0.931, 0.805, 0.55, 0.574, 0.895, 0.534, 0.742, 0.896, 0.605, 0.344, 0.99, 0.83, 0.602, 0.512, 0.773, 0.911, 0.884, 0.959, 0.918, 0.462, 0.975, 0.554, 0.969, 0.927, 0.636, 0.928, 0.742, 0.958, 0.966, 0.664, 0.956, 0.701, 0.947, 0.578, 0.594, 0.558, 0.965, 0.952, 0.974, 0.945, 0.594, 0.896, 0.976, 0.886, 0.872, 0.953, 0.894, 0.752, 0.58, 0.557, 0.983, 0.332, 0.767, 0.771, 0.383, 0.928, 0.634, 0.923, 0.982, 0.727, 0.626, 0.982, 0.982, 0.74, 0.276, 0.969, 0.947]
global origin = 1
global destination = 50