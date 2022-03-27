global arcs = [1 7; 1 41; 2 17; 2 24; 2 25; 2 30; 2 47; 3 10; 3 33; 3 37; 3 40; 3 43; 3 50; 4 9; 4 10; 4 23; 4 38; 4 45; 5 33; 5 45; 6 15; 6 23; 6 40; 6 43; 6 48; 7 4; 7 12; 7 27; 7 33; 7 34; 7 47; 8 2; 8 9; 8 12; 8 15; 8 35; 9 5; 9 12; 9 20; 9 46; 9 50; 10 13; 10 20; 10 22; 10 24; 10 33; 10 34; 10 43; 11 10; 11 13; 11 15; 11 29; 11 36; 11 43; 12 42; 12 47; 13 12; 13 34; 14 10; 14 19; 14 26; 14 28; 14 30; 14 41; 14 49; 15 5; 15 6; 15 12; 15 17; 15 29; 16 20; 16 26; 16 30; 16 34; 16 39; 16 45; 17 7; 17 11; 17 18; 17 28; 17 29; 18 3; 18 4; 18 25; 18 39; 18 43; 18 45; 18 50; 19 4; 19 5; 19 8; 19 9; 19 24; 19 25; 19 34; 19 39; 19 42; 19 44; 19 45; 19 47; 19 48; 20 11; 20 12; 20 24; 20 33; 20 50; 21 11; 21 17; 21 28; 21 45; 21 50; 22 7; 23 7; 23 13; 23 15; 23 16; 23 22; 23 36; 23 38; 23 42; 24 6; 24 17; 24 22; 24 27; 24 32; 24 38; 24 49; 25 6; 25 7; 25 14; 25 29; 25 38; 25 43; 26 24; 26 35; 26 38; 26 44; 26 45; 27 4; 27 15; 27 39; 27 44; 27 47; 28 21; 28 27; 28 35; 28 40; 28 46; 28 50; 29 23; 30 12; 30 14; 30 46; 31 24; 31 26; 31 46; 32 3; 32 18; 32 19; 32 23; 32 40; 32 46; 33 6; 33 7; 33 15; 33 25; 33 31; 33 43; 33 48; 34 6; 34 24; 34 30; 35 17; 35 20; 35 33; 35 45; 35 49; 36 11; 36 15; 36 41; 36 43; 37 9; 37 31; 37 32; 37 43; 38 8; 38 12; 38 19; 38 22; 38 32; 38 36; 38 41; 38 47; 39 3; 39 9; 39 28; 39 40; 40 3; 40 32; 41 5; 41 27; 41 36; 41 46; 42 6; 42 16; 42 36; 42 44; 42 47; 42 50; 43 8; 43 9; 43 13; 43 24; 43 48; 44 19; 44 23; 44 24; 44 37; 44 50; 45 29; 45 43; 46 4; 46 8; 46 33; 46 36; 46 39; 46 44; 47 23; 47 29; 47 30; 47 49; 48 23; 48 32; 48 34; 48 46; 48 47; 49 7; 49 14; 49 20; 49 24; 49 28; 49 34]
global d_x = [9.0, 9.0, 4.0, 8.0, 2.0, 10.0, 10.0, 1.0, 8.0, 6.0, 3.0, 4.0, 6.0, 6.0, 1.0, 4.0, 2.0, 4.0, 9.0, 10.0, 2.0, 4.0, 9.0, 1.0, 8.0, 4.0, 5.0, 2.0, 8.0, 3.0, 2.0, 4.0, 1.0, 8.0, 2.0, 4.0, 4.0, 7.0, 1.0, 4.0, 4.0, 9.0, 10.0, 3.0, 7.0, 10.0, 10.0, 4.0, 2.0, 4.0, 5.0, 6.0, 3.0, 2.0, 4.0, 4.0, 2.0, 8.0, 6.0, 1.0, 6.0, 2.0, 10.0, 4.0, 6.0, 3.0, 6.0, 9.0, 10.0, 5.0, 5.0, 6.0, 10.0, 9.0, 4.0, 6.0, 6.0, 10.0, 5.0, 9.0, 8.0, 10.0, 3.0, 6.0, 10.0, 5.0, 5.0, 10.0, 7.0, 2.0, 2.0, 10.0, 5.0, 10.0, 8.0, 1.0, 2.0, 10.0, 10.0, 10.0, 8.0, 6.0, 3.0, 4.0, 2.0, 9.0, 8.0, 4.0, 8.0, 6.0, 1.0, 9.0, 6.0, 3.0, 3.0, 2.0, 8.0, 3.0, 7.0, 7.0, 4.0, 6.0, 8.0, 10.0, 8.0, 9.0, 5.0, 4.0, 9.0, 4.0, 9.0, 3.0, 10.0, 2.0, 1.0, 8.0, 9.0, 1.0, 6.0, 1.0, 4.0, 8.0, 9.0, 4.0, 5.0, 7.0, 9.0, 5.0, 9.0, 10.0, 7.0, 5.0, 6.0, 1.0, 1.0, 5.0, 2.0, 7.0, 1.0, 5.0, 9.0, 2.0, 4.0, 2.0, 3.0, 10.0, 3.0, 3.0, 3.0, 6.0, 1.0, 6.0, 2.0, 6.0, 5.0, 2.0, 5.0, 6.0, 6.0, 5.0, 3.0, 1.0, 7.0, 8.0, 8.0, 1.0, 1.0, 1.0, 4.0, 3.0, 3.0, 5.0, 9.0, 7.0, 9.0, 1.0, 6.0, 8.0, 1.0, 8.0, 8.0, 6.0, 5.0, 10.0, 4.0, 3.0, 9.0, 4.0, 5.0, 4.0, 3.0, 4.0, 2.0, 2.0, 2.0, 6.0, 3.0, 5.0, 5.0, 6.0, 3.0, 6.0, 10.0, 10.0, 7.0, 8.0, 5.0, 5.0, 1.0, 2.0, 6.0, 5.0, 10.0, 3.0, 1.0, 8.0, 7.0, 3.0, 1.0, 5.0, 3.0, 4.0]
global b_x = 5
global d_y = [7.0, 3.0, 6.0, 6.0, 1.0, 9.0, 6.0, 10.0, 3.0, 2.0, 2.0, 8.0, 6.0, 6.0, 4.0, 2.0, 1.0, 1.0, 8.0, 10.0, 1.0, 10.0, 3.0, 7.0, 8.0, 6.0, 2.0, 2.0, 2.0, 3.0, 4.0, 8.0, 2.0, 7.0, 10.0, 8.0, 5.0, 4.0, 3.0, 10.0, 3.0, 5.0, 9.0, 6.0, 10.0, 7.0, 6.0, 9.0, 7.0, 3.0, 5.0, 2.0, 8.0, 10.0, 2.0, 4.0, 5.0, 10.0, 10.0, 3.0, 6.0, 8.0, 8.0, 3.0, 2.0, 7.0, 3.0, 2.0, 4.0, 4.0, 10.0, 2.0, 4.0, 10.0, 10.0, 3.0, 4.0, 6.0, 8.0, 5.0, 3.0, 9.0, 6.0, 1.0, 3.0, 10.0, 9.0, 2.0, 9.0, 4.0, 4.0, 6.0, 4.0, 4.0, 1.0, 1.0, 2.0, 2.0, 5.0, 4.0, 1.0, 3.0, 9.0, 7.0, 4.0, 1.0, 4.0, 1.0, 5.0, 5.0, 6.0, 8.0, 6.0, 8.0, 9.0, 6.0, 10.0, 6.0, 3.0, 1.0, 2.0, 1.0, 8.0, 6.0, 1.0, 5.0, 3.0, 8.0, 3.0, 1.0, 4.0, 5.0, 6.0, 7.0, 2.0, 7.0, 1.0, 3.0, 1.0, 6.0, 5.0, 3.0, 8.0, 6.0, 7.0, 2.0, 2.0, 7.0, 4.0, 4.0, 8.0, 3.0, 3.0, 2.0, 10.0, 6.0, 6.0, 1.0, 5.0, 8.0, 2.0, 10.0, 1.0, 4.0, 10.0, 6.0, 1.0, 2.0, 3.0, 9.0, 8.0, 4.0, 7.0, 10.0, 9.0, 7.0, 7.0, 7.0, 4.0, 7.0, 4.0, 7.0, 10.0, 2.0, 4.0, 6.0, 2.0, 2.0, 4.0, 8.0, 3.0, 2.0, 8.0, 3.0, 5.0, 10.0, 9.0, 8.0, 3.0, 5.0, 10.0, 5.0, 7.0, 9.0, 2.0, 2.0, 8.0, 5.0, 1.0, 10.0, 4.0, 9.0, 5.0, 1.0, 1.0, 7.0, 2.0, 7.0, 6.0, 2.0, 9.0, 8.0, 2.0, 5.0, 6.0, 9.0, 6.0, 1.0, 2.0, 1.0, 4.0, 7.0, 2.0, 8.0, 8.0, 5.0, 4.0, 4.0, 7.0, 2.0, 8.0, 1.0]
global b_y = 10
global p = [0.147, 0.986, 0.463, 0.911, 0.065, 0.706, 0.144, 0.819, 0.278, 0.476, 0.381, 0.315, 0.849, 0.847, 0.206, 0.717, 0.583, 0.168, 0.016, 0.744, 0.869, 0.204, 0.418, 0.423, 0.604, 0.606, 0.88, 0.467, 0.541, 0.578, 0.475, 0.452, 0.968, 0.589, 0.469, 0.002, 0.443, 0.542, 0.966, 0.365, 0.242, 0.973, 0.554, 0.279, 0.987, 0.291, 0.34, 0.167, 0.106, 0.993, 0.821, 0.059, 0.54, 0.834, 0.629, 0.255, 0.96, 0.121, 0.352, 0.9, 0.446, 0.399, 0.3, 0.57, 0.994, 0.93, 0.943, 0.822, 0.381, 0.065, 0.39, 0.572, 0.832, 0.461, 0.826, 0.593, 0.196, 0.058, 0.913, 0.744, 0.534, 0.621, 0.462, 0.51, 0.686, 0.944, 0.336, 0.663, 0.612, 0.822, 0.572, 0.441, 0.496, 0.535, 0.484, 0.643, 0.298, 0.08, 0.492, 0.755, 0.969, 0.788, 0.048, 0.171, 0.853, 0.51, 0.165, 0.033, 0.165, 0.31, 0.655, 0.958, 0.28, 0.585, 0.204, 0.743, 0.222, 0.916, 0.821, 0.383, 0.773, 0.646, 0.053, 0.169, 0.889, 0.799, 0.433, 0.348, 0.887, 0.734, 0.814, 0.755, 0.497, 0.513, 0.707, 0.626, 0.926, 0.301, 0.685, 0.336, 0.624, 0.427, 0.269, 0.038, 0.466, 0.288, 0.759, 0.607, 0.99, 0.49, 0.817, 0.052, 0.967, 0.464, 0.223, 0.455, 0.799, 0.114, 0.496, 0.099, 0.746, 0.219, 0.819, 0.731, 0.836, 0.35, 0.079, 0.401, 0.903, 0.305, 0.738, 0.027, 0.269, 0.733, 0.466, 0.909, 0.933, 0.649, 0.157, 0.967, 0.904, 0.751, 0.903, 0.503, 0.811, 0.127, 0.249, 0.75, 0.886, 0.299, 0.358, 0.737, 0.007, 0.465, 0.837, 0.118, 0.821, 0.434, 0.825, 0.78, 0.132, 0.408, 0.593, 0.308, 0.217, 0.46, 0.955, 0.368, 0.052, 0.069, 0.221, 0.549, 0.961, 0.476, 0.909, 0.83, 0.849, 0.249, 0.928, 0.747, 0.326, 0.298, 0.707, 0.986, 0.14, 0.221, 0.546, 0.833, 0.254, 0.42, 0.457, 0.249, 0.892, 0.327, 0.447, 0.779, 0.356, 0.258, 0.149, 0.665, 0.613, 0.148]
global q = [0.447, 0.99, 0.729, 0.918, 0.572, 0.942, 0.846, 0.964, 0.692, 0.834, 0.691, 0.727, 0.884, 0.912, 0.212, 0.855, 0.954, 0.914, 0.239, 0.769, 0.928, 0.812, 0.441, 0.683, 0.71, 0.91, 0.919, 0.905, 0.633, 0.633, 0.908, 0.757, 0.992, 0.883, 0.821, 0.761, 0.513, 0.925, 0.974, 0.39, 0.525, 0.981, 0.591, 0.394, 0.991, 0.74, 0.929, 0.623, 0.984, 0.993, 0.924, 0.96, 0.633, 0.916, 0.657, 0.635, 0.974, 0.709, 0.954, 0.922, 0.495, 0.437, 0.762, 0.866, 0.998, 0.964, 0.943, 0.962, 0.728, 0.613, 0.9, 0.972, 0.857, 0.742, 0.922, 0.932, 0.962, 0.112, 0.924, 0.872, 0.653, 0.7, 0.684, 0.774, 0.926, 0.975, 0.502, 0.84, 0.918, 0.824, 0.705, 0.673, 0.556, 0.965, 0.803, 0.87, 0.552, 0.963, 0.631, 0.925, 0.973, 0.8, 0.486, 0.894, 0.947, 0.611, 0.85, 0.218, 0.959, 0.553, 0.888, 0.96, 0.291, 0.628, 0.433, 0.788, 0.45, 0.94, 0.956, 0.743, 0.822, 0.827, 0.766, 0.358, 0.923, 0.929, 0.895, 0.473, 0.935, 0.932, 0.832, 0.841, 0.575, 0.855, 0.847, 0.703, 0.971, 0.733, 0.835, 0.398, 0.937, 0.475, 0.512, 0.262, 0.494, 0.801, 0.803, 0.887, 0.999, 0.53, 0.911, 0.62, 0.979, 0.994, 0.459, 0.744, 0.998, 0.872, 0.844, 0.231, 0.791, 0.809, 0.895, 0.802, 0.857, 0.855, 0.62, 0.91, 0.932, 0.816, 0.869, 0.25, 0.298, 0.76, 0.525, 0.957, 0.96, 0.721, 0.23, 0.989, 0.927, 0.904, 0.976, 0.593, 0.885, 0.982, 0.55, 0.823, 0.905, 0.925, 0.436, 0.8, 0.632, 0.974, 0.91, 0.938, 0.97, 0.949, 0.884, 0.843, 0.725, 0.513, 0.784, 0.799, 0.471, 0.917, 0.975, 0.437, 0.223, 0.601, 0.545, 0.681, 0.975, 0.894, 0.98, 0.85, 0.923, 0.843, 0.972, 0.748, 0.835, 0.852, 0.884, 0.987, 0.726, 0.499, 0.885, 0.879, 0.919, 0.971, 0.739, 0.773, 0.924, 0.394, 0.908, 0.976, 0.39, 0.887, 0.177, 0.674, 0.684, 0.521]
global origin = 1
global destination = 50