global arcs = [1 17; 1 26; 1 35; 1 41; 1 48; 2 15; 2 16; 2 17; 2 24; 2 32; 2 42; 2 45; 2 46; 3 11; 3 34; 3 36; 4 16; 4 25; 4 27; 4 31; 4 40; 4 41; 5 10; 5 14; 5 45; 6 5; 6 17; 6 23; 6 27; 6 31; 6 40; 6 50; 7 14; 7 46; 8 4; 8 9; 8 22; 8 40; 8 42; 9 14; 9 27; 9 29; 9 30; 9 37; 9 39; 9 41; 9 42; 9 48; 10 4; 10 9; 10 11; 10 21; 10 22; 10 33; 10 34; 11 2; 11 4; 11 6; 11 15; 11 30; 11 44; 12 14; 12 19; 12 32; 12 40; 12 47; 13 4; 13 5; 13 7; 13 34; 13 48; 14 3; 14 7; 14 37; 14 44; 15 9; 15 11; 15 16; 15 22; 16 29; 16 35; 16 36; 16 38; 16 49; 17 12; 17 40; 18 17; 18 29; 18 30; 18 31; 18 48; 19 8; 19 11; 19 15; 19 17; 20 3; 20 13; 20 21; 21 9; 21 15; 21 20; 21 25; 21 34; 22 34; 22 39; 23 3; 23 4; 23 11; 23 12; 23 13; 23 19; 23 22; 23 25; 23 43; 23 47; 24 20; 24 22; 24 35; 24 36; 24 37; 24 48; 25 5; 25 13; 25 18; 25 20; 25 22; 25 28; 26 10; 26 13; 26 24; 27 2; 27 14; 27 39; 27 41; 28 22; 28 35; 28 40; 29 8; 29 16; 29 22; 29 27; 29 33; 29 37; 29 41; 30 3; 30 16; 30 21; 30 38; 30 42; 31 3; 31 13; 31 16; 31 24; 31 26; 31 33; 32 4; 32 19; 32 23; 32 24; 32 28; 32 29; 32 31; 32 37; 32 50; 33 15; 33 21; 33 27; 33 41; 33 48; 34 2; 34 17; 35 11; 35 15; 35 24; 36 3; 36 5; 36 23; 36 26; 36 28; 36 37; 36 40; 37 13; 37 19; 37 24; 37 30; 37 45; 37 47; 37 49; 38 11; 38 20; 38 31; 38 37; 39 2; 39 10; 39 13; 39 21; 39 25; 39 29; 39 32; 39 33; 39 42; 40 14; 40 31; 40 39; 40 43; 41 10; 41 15; 41 20; 41 26; 41 33; 41 39; 42 9; 42 13; 42 15; 42 23; 42 37; 42 43; 42 48; 43 2; 43 4; 43 12; 43 23; 43 27; 44 3; 45 8; 45 23; 45 41; 45 43; 46 17; 46 20; 46 29; 47 13; 48 13; 48 20; 48 25; 48 27; 48 29; 48 39; 48 47; 49 12; 49 27; 49 39]
global d_x = [6.0, 6.0, 2.0, 5.0, 1.0, 5.0, 7.0, 7.0, 1.0, 6.0, 7.0, 5.0, 3.0, 5.0, 7.0, 4.0, 7.0, 2.0, 10.0, 7.0, 6.0, 10.0, 5.0, 8.0, 6.0, 1.0, 6.0, 5.0, 1.0, 2.0, 4.0, 4.0, 5.0, 9.0, 4.0, 2.0, 1.0, 2.0, 3.0, 1.0, 6.0, 9.0, 9.0, 7.0, 5.0, 2.0, 2.0, 5.0, 5.0, 5.0, 3.0, 6.0, 5.0, 5.0, 3.0, 6.0, 1.0, 6.0, 1.0, 8.0, 7.0, 8.0, 5.0, 4.0, 1.0, 5.0, 5.0, 1.0, 5.0, 3.0, 10.0, 3.0, 3.0, 7.0, 5.0, 8.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 10.0, 8.0, 8.0, 6.0, 5.0, 10.0, 5.0, 10.0, 6.0, 1.0, 10.0, 10.0, 5.0, 2.0, 4.0, 8.0, 1.0, 5.0, 3.0, 10.0, 1.0, 3.0, 3.0, 9.0, 7.0, 10.0, 6.0, 4.0, 9.0, 1.0, 2.0, 4.0, 1.0, 4.0, 1.0, 3.0, 3.0, 8.0, 9.0, 2.0, 9.0, 6.0, 7.0, 7.0, 7.0, 5.0, 5.0, 2.0, 9.0, 8.0, 9.0, 1.0, 4.0, 6.0, 5.0, 8.0, 6.0, 3.0, 9.0, 8.0, 2.0, 6.0, 6.0, 9.0, 2.0, 7.0, 8.0, 7.0, 2.0, 4.0, 9.0, 7.0, 7.0, 10.0, 7.0, 1.0, 10.0, 2.0, 1.0, 8.0, 10.0, 10.0, 2.0, 6.0, 1.0, 9.0, 10.0, 7.0, 5.0, 3.0, 4.0, 2.0, 4.0, 3.0, 1.0, 7.0, 6.0, 3.0, 9.0, 9.0, 4.0, 6.0, 10.0, 4.0, 6.0, 9.0, 9.0, 6.0, 9.0, 2.0, 5.0, 3.0, 8.0, 5.0, 3.0, 5.0, 6.0, 2.0, 9.0, 10.0, 7.0, 7.0, 6.0, 6.0, 6.0, 9.0, 1.0, 6.0, 2.0, 10.0, 4.0, 7.0, 10.0, 4.0, 7.0, 9.0, 7.0, 2.0, 1.0, 2.0, 3.0, 10.0, 7.0, 3.0, 2.0, 2.0, 9.0, 5.0, 10.0, 6.0, 3.0, 10.0, 3.0, 5.0, 6.0, 3.0, 10.0, 1.0, 10.0, 3.0]
global b_x = 5
global d_y = [1.0, 5.0, 3.0, 2.0, 5.0, 8.0, 4.0, 5.0, 9.0, 9.0, 5.0, 5.0, 6.0, 6.0, 4.0, 6.0, 8.0, 4.0, 9.0, 3.0, 8.0, 7.0, 4.0, 3.0, 1.0, 5.0, 6.0, 10.0, 2.0, 10.0, 6.0, 7.0, 6.0, 7.0, 8.0, 7.0, 6.0, 7.0, 2.0, 8.0, 4.0, 2.0, 9.0, 6.0, 8.0, 8.0, 2.0, 3.0, 5.0, 9.0, 3.0, 5.0, 10.0, 9.0, 9.0, 4.0, 7.0, 5.0, 8.0, 1.0, 7.0, 3.0, 10.0, 3.0, 9.0, 1.0, 4.0, 7.0, 1.0, 6.0, 2.0, 10.0, 5.0, 5.0, 10.0, 5.0, 9.0, 3.0, 4.0, 3.0, 8.0, 1.0, 10.0, 5.0, 5.0, 9.0, 9.0, 10.0, 5.0, 8.0, 7.0, 3.0, 8.0, 8.0, 6.0, 10.0, 9.0, 3.0, 2.0, 6.0, 4.0, 10.0, 5.0, 2.0, 3.0, 7.0, 6.0, 7.0, 1.0, 3.0, 10.0, 9.0, 4.0, 9.0, 7.0, 2.0, 2.0, 6.0, 7.0, 7.0, 4.0, 6.0, 5.0, 3.0, 7.0, 7.0, 3.0, 4.0, 6.0, 3.0, 5.0, 10.0, 10.0, 3.0, 10.0, 6.0, 4.0, 10.0, 1.0, 6.0, 8.0, 6.0, 2.0, 10.0, 6.0, 5.0, 1.0, 6.0, 10.0, 9.0, 7.0, 1.0, 3.0, 1.0, 4.0, 5.0, 2.0, 9.0, 10.0, 10.0, 9.0, 8.0, 8.0, 5.0, 4.0, 3.0, 8.0, 5.0, 10.0, 10.0, 10.0, 6.0, 2.0, 8.0, 8.0, 3.0, 2.0, 10.0, 7.0, 7.0, 8.0, 5.0, 1.0, 2.0, 10.0, 7.0, 8.0, 7.0, 3.0, 8.0, 8.0, 9.0, 4.0, 10.0, 6.0, 7.0, 2.0, 6.0, 2.0, 9.0, 3.0, 7.0, 4.0, 7.0, 7.0, 2.0, 10.0, 8.0, 2.0, 6.0, 6.0, 6.0, 5.0, 10.0, 4.0, 9.0, 4.0, 10.0, 5.0, 1.0, 5.0, 8.0, 2.0, 2.0, 5.0, 2.0, 8.0, 2.0, 4.0, 2.0, 7.0, 6.0, 2.0, 6.0, 7.0, 8.0, 6.0, 7.0, 5.0, 2.0, 1.0, 5.0]
global b_y = 10
global p = [0.022, 0.856, 0.855, 0.663, 0.405, 0.577, 0.953, 0.756, 0.087, 0.503, 0.573, 0.082, 0.205, 0.277, 0.149, 0.476, 0.524, 0.253, 0.991, 0.892, 0.699, 0.644, 0.652, 0.613, 0.746, 0.833, 0.84, 0.333, 0.489, 0.371, 0.303, 0.905, 0.088, 0.131, 0.617, 0.839, 0.572, 0.44, 0.102, 0.877, 0.835, 0.951, 0.427, 0.57, 0.746, 0.381, 0.317, 0.437, 0.887, 0.408, 0.461, 0.498, 0.964, 0.419, 0.426, 0.625, 0.955, 0.307, 0.764, 0.428, 0.525, 0.851, 0.424, 0.318, 0.057, 0.331, 0.695, 0.353, 0.09, 0.123, 0.882, 0.881, 0.06, 0.829, 0.778, 0.023, 0.059, 0.28, 0.079, 0.576, 0.097, 0.353, 0.611, 0.394, 0.478, 0.084, 0.776, 0.456, 0.828, 0.961, 0.94, 0.395, 0.155, 0.328, 0.237, 0.809, 0.411, 0.028, 0.886, 0.71, 0.241, 0.155, 0.819, 0.434, 0.542, 0.417, 0.142, 0.168, 0.512, 0.287, 0.951, 0.615, 0.156, 0.556, 0.639, 0.467, 0.022, 0.99, 0.736, 0.127, 0.49, 0.382, 0.397, 0.961, 0.924, 0.246, 0.443, 0.781, 0.29, 0.213, 0.073, 0.326, 0.925, 0.115, 0.707, 0.543, 0.136, 0.847, 0.947, 0.297, 0.466, 0.367, 0.642, 0.768, 0.206, 0.367, 0.41, 0.516, 0.927, 0.136, 0.764, 0.589, 0.06, 0.544, 0.61, 0.112, 0.219, 0.868, 0.182, 0.192, 0.139, 0.307, 0.787, 0.718, 0.599, 0.188, 0.001, 0.813, 0.669, 0.834, 0.054, 0.881, 0.194, 0.349, 0.021, 0.189, 0.406, 0.901, 0.454, 0.246, 0.025, 0.982, 0.02, 0.198, 0.942, 0.312, 0.653, 0.97, 0.09, 0.251, 0.609, 0.007, 0.287, 0.924, 0.288, 0.973, 0.432, 0.116, 0.065, 0.463, 0.556, 0.64, 0.922, 0.216, 0.466, 0.221, 0.846, 0.745, 0.86, 0.484, 0.386, 0.275, 0.038, 0.565, 0.793, 0.37, 0.697, 0.979, 0.533, 0.966, 0.902, 0.099, 0.343, 0.656, 0.236, 0.23, 0.526, 0.375, 0.5, 0.361, 0.475, 0.578, 0.409, 0.381, 0.763, 0.757, 0.759, 0.16, 0.208, 0.231, 0.759, 0.717]
global q = [0.354, 0.93, 0.937, 0.84, 0.41, 0.781, 0.993, 0.935, 0.152, 0.986, 0.885, 0.614, 0.691, 0.937, 0.443, 0.546, 0.859, 0.865, 0.993, 0.997, 0.755, 0.655, 0.831, 0.944, 0.804, 0.848, 0.852, 0.656, 0.567, 0.56, 0.322, 0.951, 0.437, 0.886, 0.932, 0.842, 0.983, 0.491, 0.135, 0.944, 0.9, 0.959, 0.596, 0.646, 0.751, 0.415, 0.354, 0.531, 0.974, 0.689, 0.96, 0.52, 0.998, 0.705, 0.681, 0.669, 0.969, 0.619, 0.914, 0.937, 0.531, 0.917, 0.532, 0.674, 0.868, 0.71, 0.737, 0.586, 0.903, 0.775, 0.965, 0.911, 0.411, 0.976, 0.979, 0.357, 0.391, 0.535, 0.489, 0.87, 0.801, 0.521, 0.655, 0.568, 0.797, 0.347, 0.994, 0.644, 0.899, 0.982, 0.992, 0.687, 0.345, 0.742, 0.809, 0.838, 0.719, 0.976, 0.988, 0.89, 0.975, 0.587, 0.885, 0.683, 0.747, 0.95, 0.671, 0.412, 0.594, 0.99, 0.994, 0.705, 0.729, 0.792, 0.695, 0.816, 0.195, 0.993, 0.953, 0.741, 0.555, 0.859, 0.919, 0.962, 0.93, 0.525, 0.85, 0.979, 0.291, 0.837, 0.138, 0.909, 0.95, 0.861, 0.953, 0.803, 0.921, 0.887, 0.988, 0.493, 0.567, 0.935, 0.972, 0.953, 0.938, 0.909, 0.819, 0.953, 0.956, 0.139, 0.818, 0.911, 0.959, 0.844, 0.62, 0.356, 0.618, 0.886, 0.683, 0.608, 0.453, 0.882, 0.9, 0.736, 0.93, 0.394, 0.273, 0.846, 0.689, 0.948, 0.992, 0.917, 0.453, 0.446, 0.091, 0.418, 0.486, 0.932, 0.617, 0.96, 0.388, 0.997, 0.193, 0.351, 0.967, 0.541, 0.8, 0.997, 0.272, 0.426, 0.745, 0.188, 0.755, 0.954, 0.479, 0.993, 0.604, 0.889, 0.574, 0.608, 0.904, 0.99, 0.998, 0.408, 0.953, 0.957, 0.927, 0.859, 0.966, 0.788, 0.731, 0.816, 0.766, 0.954, 0.955, 0.529, 0.779, 0.999, 0.65, 0.98, 0.929, 0.866, 0.349, 0.745, 0.505, 0.996, 0.755, 0.692, 0.512, 0.531, 0.601, 0.647, 0.492, 0.382, 0.816, 0.9, 0.834, 0.957, 0.821, 0.584, 0.964, 0.948]
global origin = 1
global destination = 50