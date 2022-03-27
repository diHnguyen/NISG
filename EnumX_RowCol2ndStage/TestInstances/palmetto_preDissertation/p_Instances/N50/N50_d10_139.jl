global arcs = [1 8; 1 30; 1 35; 1 39; 1 48; 2 7; 2 8; 2 13; 2 21; 2 24; 2 26; 2 32; 2 38; 2 40; 3 2; 3 9; 3 33; 3 34; 3 42; 4 9; 4 19; 4 20; 4 31; 4 39; 4 49; 5 6; 5 13; 5 22; 5 26; 5 39; 5 40; 6 4; 6 22; 6 23; 6 25; 6 27; 6 39; 6 40; 7 2; 7 13; 7 25; 7 42; 8 5; 8 13; 8 24; 8 29; 8 30; 9 2; 9 4; 9 12; 10 8; 10 14; 10 16; 10 24; 10 31; 10 32; 10 33; 10 38; 10 47; 10 48; 10 50; 11 7; 11 17; 11 34; 11 37; 11 48; 12 11; 12 15; 12 26; 12 28; 12 31; 12 32; 12 38; 12 40; 12 46; 12 47; 12 48; 12 50; 13 19; 13 26; 13 30; 13 38; 13 41; 13 42; 13 48; 14 4; 14 39; 14 42; 15 8; 15 11; 15 28; 15 42; 16 8; 16 9; 16 18; 16 22; 16 25; 17 25; 17 30; 18 4; 18 23; 19 5; 19 16; 19 25; 19 26; 19 30; 20 21; 20 32; 20 43; 20 47; 21 14; 21 20; 21 42; 22 3; 22 12; 22 24; 22 27; 22 48; 23 12; 23 22; 23 32; 23 36; 23 38; 23 39; 23 43; 24 2; 24 10; 24 11; 24 12; 24 14; 24 18; 24 42; 24 48; 25 11; 25 28; 26 8; 26 10; 26 36; 26 44; 27 8; 27 13; 27 26; 27 41; 27 47; 28 6; 29 4; 30 6; 30 16; 30 21; 30 25; 30 32; 30 37; 30 38; 30 44; 31 8; 31 15; 31 36; 32 13; 32 18; 32 29; 32 45; 33 25; 33 28; 33 39; 34 2; 34 7; 34 23; 34 24; 34 29; 34 49; 35 3; 35 4; 35 21; 35 28; 35 39; 35 49; 36 8; 36 13; 36 16; 36 17; 36 35; 36 42; 37 21; 38 3; 38 5; 38 16; 38 26; 38 48; 39 13; 39 16; 40 45; 40 48; 41 10; 42 8; 42 14; 42 19; 42 21; 42 32; 42 45; 43 8; 43 9; 43 10; 43 22; 43 31; 43 37; 44 10; 44 14; 44 21; 44 27; 44 34; 44 46; 45 31; 46 2; 46 25; 47 2; 47 11; 47 17; 48 17; 48 20; 48 25; 48 27; 48 35; 48 39; 49 4; 49 17; 49 41]
global d_x = [4.0, 1.0, 4.0, 1.0, 7.0, 8.0, 5.0, 3.0, 4.0, 1.0, 10.0, 3.0, 4.0, 3.0, 8.0, 7.0, 1.0, 3.0, 3.0, 1.0, 1.0, 2.0, 3.0, 4.0, 3.0, 4.0, 6.0, 8.0, 7.0, 4.0, 1.0, 10.0, 9.0, 1.0, 4.0, 1.0, 2.0, 5.0, 3.0, 6.0, 1.0, 7.0, 9.0, 4.0, 10.0, 2.0, 8.0, 4.0, 2.0, 10.0, 2.0, 10.0, 8.0, 8.0, 9.0, 8.0, 1.0, 8.0, 9.0, 9.0, 9.0, 2.0, 6.0, 8.0, 9.0, 5.0, 7.0, 1.0, 6.0, 7.0, 1.0, 6.0, 1.0, 3.0, 6.0, 1.0, 8.0, 5.0, 9.0, 9.0, 10.0, 3.0, 9.0, 4.0, 10.0, 7.0, 4.0, 4.0, 2.0, 7.0, 5.0, 3.0, 3.0, 9.0, 6.0, 2.0, 8.0, 3.0, 5.0, 8.0, 6.0, 1.0, 6.0, 4.0, 3.0, 6.0, 7.0, 4.0, 7.0, 4.0, 7.0, 8.0, 1.0, 4.0, 3.0, 7.0, 1.0, 9.0, 5.0, 5.0, 4.0, 9.0, 2.0, 5.0, 5.0, 5.0, 10.0, 1.0, 2.0, 4.0, 6.0, 5.0, 8.0, 8.0, 8.0, 2.0, 5.0, 8.0, 9.0, 7.0, 4.0, 3.0, 1.0, 4.0, 1.0, 7.0, 3.0, 9.0, 10.0, 4.0, 6.0, 1.0, 4.0, 8.0, 9.0, 1.0, 6.0, 5.0, 7.0, 2.0, 4.0, 3.0, 7.0, 1.0, 6.0, 1.0, 3.0, 4.0, 3.0, 1.0, 5.0, 6.0, 2.0, 7.0, 1.0, 2.0, 7.0, 8.0, 1.0, 6.0, 9.0, 8.0, 2.0, 4.0, 5.0, 5.0, 10.0, 3.0, 7.0, 3.0, 3.0, 3.0, 8.0, 2.0, 8.0, 10.0, 3.0, 2.0, 5.0, 3.0, 8.0, 5.0, 10.0, 10.0, 4.0, 4.0, 2.0, 7.0, 8.0, 4.0, 3.0, 5.0, 6.0, 2.0, 2.0, 9.0, 4.0, 2.0, 7.0, 6.0, 1.0, 2.0, 7.0, 7.0, 4.0, 10.0]
global b_x = 5
global d_y = [8.0, 4.0, 9.0, 5.0, 1.0, 2.0, 9.0, 6.0, 7.0, 4.0, 1.0, 2.0, 5.0, 7.0, 4.0, 1.0, 8.0, 8.0, 2.0, 1.0, 3.0, 7.0, 2.0, 7.0, 5.0, 2.0, 8.0, 2.0, 5.0, 5.0, 7.0, 7.0, 4.0, 10.0, 6.0, 9.0, 6.0, 1.0, 8.0, 2.0, 10.0, 4.0, 4.0, 3.0, 6.0, 9.0, 5.0, 1.0, 7.0, 5.0, 5.0, 2.0, 9.0, 4.0, 2.0, 8.0, 6.0, 4.0, 2.0, 2.0, 5.0, 2.0, 4.0, 5.0, 7.0, 9.0, 7.0, 6.0, 3.0, 4.0, 7.0, 7.0, 8.0, 10.0, 1.0, 8.0, 5.0, 5.0, 2.0, 6.0, 6.0, 10.0, 8.0, 6.0, 3.0, 7.0, 3.0, 6.0, 4.0, 9.0, 1.0, 9.0, 5.0, 9.0, 6.0, 7.0, 2.0, 6.0, 4.0, 3.0, 2.0, 10.0, 10.0, 1.0, 6.0, 8.0, 1.0, 10.0, 7.0, 9.0, 6.0, 5.0, 4.0, 4.0, 7.0, 8.0, 10.0, 5.0, 4.0, 1.0, 9.0, 8.0, 4.0, 3.0, 5.0, 6.0, 8.0, 7.0, 5.0, 5.0, 2.0, 8.0, 8.0, 6.0, 9.0, 2.0, 4.0, 1.0, 3.0, 2.0, 7.0, 10.0, 2.0, 9.0, 8.0, 5.0, 2.0, 8.0, 6.0, 9.0, 9.0, 8.0, 4.0, 1.0, 1.0, 2.0, 6.0, 2.0, 5.0, 9.0, 3.0, 2.0, 4.0, 10.0, 6.0, 3.0, 7.0, 10.0, 10.0, 9.0, 7.0, 4.0, 1.0, 5.0, 1.0, 3.0, 7.0, 7.0, 3.0, 5.0, 4.0, 10.0, 10.0, 9.0, 2.0, 7.0, 3.0, 9.0, 5.0, 5.0, 5.0, 7.0, 5.0, 2.0, 7.0, 4.0, 3.0, 9.0, 8.0, 2.0, 5.0, 3.0, 4.0, 7.0, 8.0, 6.0, 4.0, 9.0, 3.0, 3.0, 1.0, 7.0, 10.0, 1.0, 2.0, 1.0, 10.0, 10.0, 6.0, 1.0, 2.0, 10.0, 5.0, 7.0, 1.0, 1.0]
global b_y = 10
global p = [0.386, 0.636, 0.644, 0.932, 0.474, 0.38, 0.494, 0.988, 0.474, 0.037, 0.04, 0.001, 0.079, 0.719, 0.335, 0.119, 0.595, 0.786, 0.014, 0.453, 0.892, 0.78, 0.257, 0.905, 0.266, 0.623, 0.112, 0.14, 0.61, 0.638, 0.094, 0.194, 0.311, 0.804, 0.476, 0.659, 0.264, 0.742, 0.141, 0.788, 0.903, 0.668, 0.914, 0.619, 0.055, 0.409, 0.832, 0.537, 0.017, 0.002, 0.584, 0.511, 0.232, 0.311, 0.117, 0.351, 0.152, 0.275, 0.537, 0.326, 0.362, 0.768, 0.713, 0.323, 0.899, 0.396, 0.971, 0.557, 0.823, 0.293, 0.497, 0.118, 0.464, 0.094, 0.097, 0.93, 0.063, 0.475, 0.587, 0.108, 0.242, 0.863, 0.901, 0.834, 0.198, 0.45, 0.04, 0.664, 0.888, 0.14, 0.579, 0.766, 0.215, 0.863, 0.245, 0.96, 0.017, 0.511, 0.628, 0.496, 0.652, 0.71, 0.367, 0.243, 0.498, 0.74, 0.009, 0.246, 0.841, 0.955, 0.385, 0.173, 0.544, 0.903, 0.607, 0.789, 0.16, 0.281, 0.877, 0.891, 0.148, 0.541, 0.247, 0.806, 0.99, 0.482, 0.497, 0.169, 0.133, 0.723, 0.024, 0.826, 0.907, 0.025, 0.493, 0.93, 0.841, 0.772, 0.189, 0.563, 0.065, 0.825, 0.883, 0.424, 0.874, 0.149, 0.798, 0.769, 0.159, 0.476, 0.686, 0.207, 0.424, 0.465, 0.468, 0.8, 0.158, 0.028, 0.302, 0.607, 0.771, 0.072, 0.275, 0.551, 0.905, 0.597, 0.176, 0.889, 0.259, 0.162, 0.443, 0.145, 0.911, 0.663, 0.109, 0.79, 0.509, 0.642, 0.647, 0.773, 0.792, 0.372, 0.932, 0.697, 0.781, 0.314, 0.074, 0.353, 0.404, 0.139, 0.412, 0.784, 0.749, 0.482, 0.397, 0.098, 0.325, 0.513, 0.837, 0.28, 0.485, 0.562, 0.123, 0.877, 0.062, 0.467, 0.242, 0.821, 0.862, 0.49, 0.721, 0.828, 0.815, 0.694, 0.685, 0.983, 0.532, 0.161, 0.009, 0.458, 0.312, 0.789, 0.431, 0.155, 0.414, 0.798]
global q = [0.873, 0.674, 0.933, 0.944, 0.892, 0.765, 0.547, 0.991, 0.474, 0.054, 0.276, 0.167, 0.23, 0.749, 0.889, 0.852, 0.901, 0.925, 0.052, 0.911, 0.953, 0.946, 0.861, 0.972, 0.871, 0.672, 0.235, 0.805, 0.923, 0.9, 0.971, 0.515, 0.4, 0.969, 0.757, 0.728, 0.75, 0.977, 0.543, 0.824, 0.984, 0.777, 0.943, 0.798, 0.677, 0.594, 0.836, 0.813, 0.596, 0.365, 0.631, 0.959, 0.816, 0.505, 0.893, 0.381, 0.834, 0.51, 0.674, 0.428, 0.791, 0.923, 0.946, 0.924, 0.979, 0.419, 0.987, 0.697, 0.935, 0.593, 0.713, 0.724, 0.864, 0.573, 0.274, 0.99, 0.135, 0.894, 0.629, 0.8, 0.855, 0.968, 0.956, 0.976, 0.449, 0.836, 0.603, 0.749, 0.968, 0.426, 0.954, 0.839, 0.248, 0.887, 0.825, 0.963, 0.073, 0.767, 0.727, 0.656, 0.929, 0.939, 0.761, 0.902, 0.522, 0.838, 0.315, 0.872, 0.977, 0.975, 0.505, 0.325, 0.853, 0.925, 0.842, 0.949, 0.703, 0.627, 0.915, 0.904, 0.787, 0.784, 0.446, 0.807, 0.99, 0.939, 0.911, 0.355, 0.96, 0.82, 0.457, 0.847, 0.982, 0.042, 0.577, 0.972, 0.951, 0.813, 0.485, 0.652, 0.496, 0.995, 0.885, 0.69, 0.907, 0.389, 0.896, 0.876, 0.169, 0.611, 0.997, 0.612, 0.72, 0.544, 0.703, 0.924, 0.561, 0.992, 0.968, 0.918, 0.831, 0.404, 0.628, 0.986, 0.973, 0.652, 0.301, 0.928, 0.387, 0.586, 0.951, 0.668, 0.923, 0.845, 0.392, 0.923, 0.673, 0.941, 0.954, 0.889, 0.959, 0.382, 0.937, 0.925, 0.984, 0.41, 0.862, 0.483, 0.645, 0.613, 0.57, 0.836, 0.967, 0.909, 0.7, 0.831, 0.972, 0.896, 0.945, 0.418, 0.705, 0.634, 0.517, 0.965, 0.826, 0.95, 0.972, 0.829, 0.925, 0.844, 0.894, 0.92, 0.948, 0.887, 0.739, 0.986, 0.815, 0.386, 0.601, 0.898, 0.794, 0.901, 0.51, 0.81, 0.643, 0.954]
global origin = 1
global destination = 50