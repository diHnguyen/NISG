global arcs = [1 20; 1 22; 1 28; 1 42; 1 44; 1 47; 2 3; 2 8; 2 41; 2 50; 3 12; 3 17; 3 31; 4 12; 4 46; 4 47; 5 4; 5 15; 5 21; 5 30; 5 45; 6 4; 6 18; 6 21; 6 24; 6 36; 6 44; 6 46; 7 6; 7 18; 8 12; 8 22; 8 33; 9 10; 9 24; 9 25; 9 26; 9 28; 10 11; 10 13; 10 34; 11 5; 11 15; 12 10; 12 21; 13 2; 13 22; 13 23; 13 26; 14 4; 14 11; 14 45; 15 38; 16 10; 16 20; 16 21; 16 23; 16 24; 16 33; 16 44; 17 29; 18 8; 18 16; 18 19; 18 20; 18 21; 18 26; 18 27; 18 33; 18 37; 18 42; 19 6; 19 24; 19 35; 20 16; 20 37; 20 42; 20 44; 21 6; 21 11; 21 20; 21 23; 21 25; 21 39; 22 2; 22 16; 22 34; 22 46; 23 13; 23 18; 23 25; 23 39; 23 41; 24 6; 24 9; 24 12; 24 13; 24 23; 24 27; 24 40; 24 44; 24 49; 25 8; 25 10; 25 11; 25 19; 26 25; 26 27; 26 49; 27 40; 28 11; 28 19; 28 33; 28 44; 28 45; 29 20; 29 21; 29 26; 29 30; 29 48; 30 6; 30 9; 30 16; 30 21; 30 23; 30 31; 30 37; 30 46; 30 49; 31 25; 31 28; 31 43; 32 14; 32 35; 33 4; 33 32; 33 36; 33 48; 34 11; 34 14; 34 19; 34 23; 34 33; 34 38; 34 49; 35 4; 35 12; 35 28; 35 31; 35 34; 35 41; 36 18; 36 43; 36 45; 36 50; 37 2; 37 8; 37 11; 37 45; 38 37; 38 47; 39 7; 39 9; 39 15; 39 27; 40 2; 40 12; 40 30; 40 31; 40 34; 40 41; 40 43; 40 49; 41 23; 41 48; 42 18; 42 21; 42 44; 43 6; 43 15; 43 23; 43 35; 43 38; 43 44; 44 10; 44 16; 44 18; 44 20; 44 33; 44 37; 44 49; 45 2; 45 12; 45 23; 45 24; 45 25; 45 50; 46 20; 46 22; 46 44; 46 45; 46 47; 47 13; 47 18; 47 20; 47 28; 47 30; 47 31; 47 42; 48 11; 49 15; 49 18; 49 20; 49 36; 49 38; 49 46; 49 50]
global d_x = [2.0, 7.0, 5.0, 4.0, 9.0, 5.0, 5.0, 4.0, 6.0, 8.0, 3.0, 6.0, 8.0, 1.0, 6.0, 9.0, 5.0, 4.0, 8.0, 6.0, 9.0, 5.0, 4.0, 5.0, 5.0, 1.0, 8.0, 9.0, 4.0, 5.0, 2.0, 2.0, 9.0, 5.0, 5.0, 1.0, 3.0, 7.0, 7.0, 10.0, 8.0, 3.0, 3.0, 4.0, 4.0, 9.0, 4.0, 7.0, 4.0, 9.0, 4.0, 2.0, 5.0, 9.0, 6.0, 7.0, 1.0, 2.0, 6.0, 5.0, 9.0, 2.0, 3.0, 9.0, 6.0, 9.0, 2.0, 1.0, 5.0, 7.0, 5.0, 1.0, 10.0, 7.0, 4.0, 4.0, 4.0, 7.0, 1.0, 5.0, 5.0, 1.0, 9.0, 1.0, 6.0, 8.0, 10.0, 3.0, 5.0, 2.0, 4.0, 4.0, 2.0, 9.0, 1.0, 9.0, 7.0, 6.0, 9.0, 4.0, 9.0, 7.0, 10.0, 1.0, 3.0, 10.0, 10.0, 3.0, 10.0, 7.0, 1.0, 5.0, 7.0, 10.0, 6.0, 4.0, 10.0, 2.0, 2.0, 8.0, 5.0, 3.0, 4.0, 6.0, 3.0, 4.0, 9.0, 7.0, 6.0, 2.0, 4.0, 9.0, 1.0, 7.0, 8.0, 3.0, 1.0, 7.0, 9.0, 7.0, 1.0, 10.0, 7.0, 1.0, 7.0, 10.0, 1.0, 2.0, 8.0, 5.0, 7.0, 1.0, 7.0, 3.0, 3.0, 3.0, 8.0, 6.0, 6.0, 8.0, 3.0, 3.0, 8.0, 5.0, 8.0, 8.0, 7.0, 9.0, 1.0, 7.0, 7.0, 3.0, 6.0, 1.0, 5.0, 6.0, 8.0, 1.0, 7.0, 7.0, 4.0, 1.0, 3.0, 4.0, 4.0, 1.0, 3.0, 5.0, 3.0, 7.0, 9.0, 2.0, 5.0, 5.0, 6.0, 8.0, 5.0, 8.0, 6.0, 9.0, 6.0, 1.0, 10.0, 3.0, 1.0, 1.0, 10.0, 4.0, 3.0, 1.0, 9.0, 7.0, 10.0, 3.0, 10.0, 9.0, 2.0]
global b_x = 5
global d_y = [10.0, 5.0, 6.0, 5.0, 5.0, 7.0, 8.0, 9.0, 3.0, 3.0, 1.0, 7.0, 5.0, 8.0, 5.0, 9.0, 10.0, 1.0, 7.0, 3.0, 5.0, 10.0, 10.0, 8.0, 10.0, 3.0, 10.0, 10.0, 5.0, 5.0, 1.0, 9.0, 9.0, 5.0, 8.0, 9.0, 2.0, 8.0, 1.0, 4.0, 8.0, 5.0, 4.0, 4.0, 1.0, 9.0, 9.0, 6.0, 3.0, 3.0, 1.0, 4.0, 1.0, 4.0, 1.0, 4.0, 9.0, 3.0, 2.0, 7.0, 7.0, 8.0, 6.0, 1.0, 2.0, 8.0, 10.0, 6.0, 10.0, 3.0, 10.0, 3.0, 5.0, 6.0, 4.0, 5.0, 1.0, 3.0, 4.0, 7.0, 7.0, 6.0, 4.0, 2.0, 10.0, 2.0, 5.0, 9.0, 6.0, 7.0, 4.0, 10.0, 9.0, 9.0, 10.0, 9.0, 7.0, 5.0, 6.0, 10.0, 6.0, 4.0, 6.0, 5.0, 3.0, 6.0, 9.0, 8.0, 1.0, 8.0, 5.0, 9.0, 8.0, 7.0, 3.0, 6.0, 10.0, 9.0, 2.0, 2.0, 4.0, 3.0, 7.0, 2.0, 8.0, 4.0, 7.0, 5.0, 10.0, 7.0, 5.0, 6.0, 7.0, 6.0, 3.0, 9.0, 8.0, 1.0, 5.0, 4.0, 2.0, 9.0, 7.0, 2.0, 7.0, 10.0, 7.0, 10.0, 9.0, 9.0, 7.0, 2.0, 7.0, 1.0, 4.0, 10.0, 5.0, 3.0, 3.0, 2.0, 7.0, 6.0, 1.0, 1.0, 1.0, 4.0, 9.0, 10.0, 2.0, 3.0, 1.0, 10.0, 9.0, 4.0, 10.0, 10.0, 7.0, 2.0, 6.0, 2.0, 6.0, 4.0, 7.0, 6.0, 3.0, 10.0, 1.0, 8.0, 1.0, 2.0, 3.0, 2.0, 5.0, 8.0, 7.0, 3.0, 9.0, 10.0, 7.0, 5.0, 5.0, 2.0, 1.0, 8.0, 8.0, 5.0, 2.0, 2.0, 6.0, 2.0, 7.0, 7.0, 7.0, 8.0, 3.0, 2.0, 8.0]
global b_y = 10
global p = [0.442, 0.144, 0.419, 0.612, 0.145, 0.225, 0.691, 0.634, 0.391, 0.001, 0.022, 0.711, 0.353, 0.28, 0.426, 0.487, 0.615, 0.367, 0.363, 0.663, 0.84, 0.189, 0.597, 0.204, 0.208, 0.04, 0.359, 0.764, 0.559, 0.962, 0.866, 0.124, 0.292, 0.092, 0.119, 0.22, 0.295, 0.156, 0.794, 0.589, 0.097, 0.595, 0.391, 0.014, 0.369, 0.816, 0.825, 0.852, 0.546, 0.228, 0.923, 0.617, 0.653, 0.432, 0.693, 0.9, 0.418, 0.873, 0.564, 0.881, 0.393, 0.454, 0.272, 0.366, 0.414, 0.08, 0.356, 0.016, 0.203, 0.914, 0.43, 0.145, 0.513, 0.79, 0.727, 0.377, 0.191, 0.104, 0.61, 0.297, 0.749, 0.509, 0.371, 0.844, 0.201, 0.918, 0.048, 0.797, 0.991, 0.264, 0.71, 0.974, 0.326, 0.862, 0.391, 0.292, 0.644, 0.111, 0.295, 0.055, 0.038, 0.839, 0.074, 0.178, 0.421, 0.829, 0.926, 0.021, 0.777, 0.435, 0.22, 0.321, 0.299, 0.241, 0.856, 0.039, 0.301, 0.958, 0.146, 0.921, 0.393, 0.546, 0.92, 0.736, 0.995, 0.177, 0.703, 0.675, 0.695, 0.896, 0.23, 0.561, 0.28, 0.979, 0.332, 0.948, 0.145, 0.386, 0.414, 0.458, 0.854, 0.582, 0.353, 0.618, 0.887, 0.568, 0.701, 0.931, 0.714, 0.201, 0.914, 0.892, 0.61, 0.232, 0.999, 0.034, 0.407, 0.692, 0.984, 0.551, 0.446, 0.083, 0.203, 0.665, 0.915, 0.818, 0.351, 0.195, 0.317, 0.188, 0.696, 0.967, 0.467, 0.366, 0.7, 0.718, 0.269, 0.574, 0.591, 0.425, 0.556, 0.183, 0.372, 0.563, 0.698, 0.665, 0.96, 0.132, 0.312, 0.716, 0.052, 0.831, 0.082, 0.498, 0.596, 0.512, 0.351, 0.299, 0.412, 0.343, 0.935, 0.257, 0.513, 0.728, 0.752, 0.603, 0.894, 0.709, 0.192, 0.511, 0.347, 0.811, 0.422, 0.041, 0.731, 0.109, 0.716]
global q = [0.832, 0.333, 0.966, 0.928, 0.201, 0.898, 0.994, 0.922, 0.439, 0.605, 0.742, 0.897, 0.95, 0.655, 0.972, 0.857, 0.793, 0.799, 0.499, 0.915, 0.98, 0.887, 0.841, 0.516, 0.923, 0.696, 0.406, 0.94, 0.919, 0.968, 0.947, 0.154, 0.662, 0.477, 0.28, 0.233, 0.824, 0.656, 0.795, 0.682, 0.668, 0.82, 0.834, 0.858, 0.773, 0.864, 0.842, 0.967, 0.98, 0.722, 0.935, 0.847, 0.674, 0.531, 0.763, 0.956, 0.682, 0.913, 0.757, 0.972, 0.694, 0.799, 0.492, 0.671, 0.813, 0.711, 0.991, 0.952, 0.804, 0.933, 0.881, 0.324, 0.779, 0.913, 0.996, 0.998, 0.288, 0.703, 0.933, 0.383, 0.831, 0.821, 0.5, 0.902, 0.41, 0.97, 0.774, 0.87, 0.992, 0.54, 0.987, 0.99, 0.983, 0.926, 0.642, 0.483, 0.883, 0.804, 0.949, 0.397, 0.806, 0.849, 0.935, 0.759, 0.56, 0.907, 0.951, 0.333, 0.918, 0.891, 0.709, 0.621, 0.538, 0.871, 0.897, 0.528, 0.706, 0.985, 0.539, 0.94, 0.606, 0.631, 0.992, 0.82, 0.999, 0.401, 0.888, 0.932, 0.884, 0.903, 0.468, 0.947, 0.665, 0.987, 0.579, 0.954, 0.579, 0.925, 0.83, 0.572, 0.862, 0.911, 0.685, 0.792, 0.965, 0.783, 0.901, 0.978, 0.76, 0.206, 0.954, 0.982, 0.718, 0.842, 0.999, 0.147, 0.719, 0.869, 0.999, 0.73, 0.623, 0.328, 0.851, 0.812, 0.97, 0.88, 0.699, 0.921, 0.492, 0.4, 0.934, 0.968, 0.749, 0.838, 0.781, 0.873, 0.957, 0.606, 0.803, 0.695, 0.67, 0.527, 0.561, 0.779, 0.737, 0.757, 0.995, 0.333, 0.833, 0.88, 0.368, 0.929, 0.322, 0.821, 0.626, 0.841, 0.708, 0.844, 0.967, 0.502, 0.958, 0.931, 0.934, 0.875, 0.955, 0.945, 0.976, 0.719, 0.888, 0.866, 0.551, 0.973, 0.902, 0.642, 0.944, 0.349, 0.734]
global origin = 1
global destination = 50