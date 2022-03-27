global arcs = [1 6; 1 15; 1 19; 1 22; 1 27; 1 32; 1 39; 2 7; 2 27; 2 33; 2 42; 2 43; 3 5; 3 26; 3 29; 3 48; 3 50; 4 16; 4 31; 4 33; 4 50; 5 19; 5 43; 5 45; 6 7; 6 13; 6 37; 6 38; 7 5; 7 10; 7 16; 7 39; 7 43; 7 49; 8 5; 8 12; 8 18; 8 33; 9 4; 9 15; 9 40; 9 43; 10 40; 10 44; 10 49; 11 3; 11 32; 11 35; 11 39; 11 40; 11 42; 12 11; 12 20; 12 24; 12 30; 12 34; 12 37; 13 11; 13 15; 13 16; 13 22; 13 26; 13 27; 13 32; 13 33; 13 40; 13 46; 14 12; 15 7; 15 12; 15 20; 15 33; 16 5; 16 26; 16 28; 16 46; 17 27; 18 3; 18 21; 18 22; 18 33; 18 50; 19 2; 19 21; 19 26; 19 29; 19 44; 19 46; 19 49; 20 5; 20 12; 20 28; 21 4; 21 10; 21 13; 21 24; 21 26; 21 41; 21 44; 21 46; 22 4; 22 8; 22 16; 22 17; 22 29; 22 37; 22 42; 22 46; 22 47; 23 5; 23 6; 23 8; 23 21; 24 3; 24 6; 24 16; 24 18; 24 36; 24 47; 24 48; 24 49; 25 7; 25 20; 25 28; 25 36; 25 46; 25 49; 26 4; 26 27; 26 38; 26 40; 26 41; 26 50; 27 5; 27 20; 27 33; 27 36; 28 14; 28 49; 29 4; 29 12; 29 17; 29 24; 29 42; 29 47; 30 4; 30 12; 30 16; 30 17; 30 27; 30 28; 30 49; 31 8; 31 14; 31 20; 31 32; 31 33; 31 36; 31 46; 32 3; 33 4; 33 14; 33 42; 34 11; 34 35; 34 47; 35 3; 35 15; 35 18; 35 45; 35 47; 36 6; 36 11; 36 17; 36 19; 36 37; 37 17; 37 23; 37 33; 37 40; 38 6; 39 2; 39 5; 39 16; 39 31; 40 20; 40 25; 40 33; 40 42; 40 46; 41 2; 41 11; 41 15; 41 16; 41 25; 42 8; 42 33; 42 44; 42 45; 43 6; 43 10; 43 13; 43 17; 43 35; 43 39; 44 7; 44 22; 44 27; 45 3; 45 16; 45 19; 45 21; 45 28; 45 42; 46 7; 46 9; 46 20; 46 25; 46 38; 46 39; 46 47; 47 8; 47 9; 47 14; 47 15; 47 19; 47 29; 47 34; 47 35; 47 41; 48 46; 49 3; 49 10; 49 31; 49 42; 49 43; 49 44; 49 50]
global d_x = [6.0, 4.0, 5.0, 6.0, 3.0, 6.0, 1.0, 6.0, 2.0, 4.0, 5.0, 5.0, 3.0, 9.0, 8.0, 3.0, 4.0, 1.0, 5.0, 7.0, 9.0, 8.0, 4.0, 5.0, 5.0, 6.0, 3.0, 10.0, 10.0, 6.0, 4.0, 2.0, 10.0, 9.0, 8.0, 10.0, 10.0, 1.0, 2.0, 9.0, 4.0, 7.0, 2.0, 1.0, 1.0, 3.0, 3.0, 5.0, 3.0, 6.0, 9.0, 4.0, 4.0, 10.0, 8.0, 4.0, 8.0, 4.0, 7.0, 8.0, 8.0, 9.0, 5.0, 7.0, 8.0, 10.0, 6.0, 10.0, 8.0, 9.0, 7.0, 2.0, 5.0, 6.0, 1.0, 2.0, 2.0, 10.0, 2.0, 5.0, 4.0, 10.0, 3.0, 5.0, 1.0, 2.0, 7.0, 9.0, 7.0, 2.0, 7.0, 8.0, 2.0, 7.0, 2.0, 10.0, 3.0, 3.0, 10.0, 4.0, 9.0, 6.0, 2.0, 6.0, 3.0, 2.0, 2.0, 5.0, 3.0, 3.0, 5.0, 9.0, 3.0, 8.0, 2.0, 7.0, 2.0, 1.0, 10.0, 3.0, 5.0, 7.0, 5.0, 3.0, 2.0, 9.0, 7.0, 7.0, 2.0, 6.0, 10.0, 6.0, 9.0, 8.0, 1.0, 9.0, 3.0, 5.0, 3.0, 10.0, 3.0, 4.0, 4.0, 10.0, 6.0, 8.0, 3.0, 6.0, 7.0, 7.0, 3.0, 9.0, 4.0, 7.0, 3.0, 6.0, 9.0, 5.0, 7.0, 6.0, 8.0, 9.0, 5.0, 9.0, 8.0, 6.0, 4.0, 5.0, 4.0, 4.0, 3.0, 10.0, 3.0, 8.0, 2.0, 1.0, 3.0, 6.0, 6.0, 5.0, 1.0, 7.0, 1.0, 8.0, 10.0, 5.0, 2.0, 9.0, 6.0, 1.0, 2.0, 2.0, 6.0, 9.0, 6.0, 7.0, 10.0, 6.0, 7.0, 4.0, 9.0, 6.0, 5.0, 10.0, 10.0, 8.0, 4.0, 4.0, 5.0, 3.0, 10.0, 2.0, 8.0, 2.0, 6.0, 1.0, 4.0, 3.0, 7.0, 7.0, 4.0, 1.0, 2.0, 6.0, 7.0, 6.0, 9.0, 8.0, 5.0, 9.0, 6.0, 1.0, 6.0, 2.0, 6.0, 9.0, 3.0, 6.0]
global b_x = 5
global d_y = [10.0, 6.0, 3.0, 9.0, 10.0, 4.0, 9.0, 4.0, 6.0, 10.0, 9.0, 7.0, 8.0, 9.0, 4.0, 6.0, 3.0, 4.0, 5.0, 1.0, 6.0, 4.0, 9.0, 5.0, 2.0, 2.0, 1.0, 3.0, 6.0, 4.0, 7.0, 3.0, 4.0, 1.0, 9.0, 6.0, 1.0, 5.0, 1.0, 4.0, 6.0, 8.0, 8.0, 7.0, 3.0, 4.0, 8.0, 8.0, 3.0, 8.0, 5.0, 2.0, 1.0, 5.0, 2.0, 2.0, 2.0, 7.0, 6.0, 1.0, 6.0, 2.0, 4.0, 9.0, 7.0, 3.0, 4.0, 3.0, 6.0, 9.0, 2.0, 7.0, 10.0, 9.0, 9.0, 8.0, 7.0, 9.0, 4.0, 4.0, 2.0, 8.0, 8.0, 8.0, 3.0, 4.0, 4.0, 2.0, 6.0, 4.0, 8.0, 5.0, 3.0, 10.0, 4.0, 4.0, 5.0, 5.0, 10.0, 10.0, 9.0, 7.0, 1.0, 9.0, 9.0, 4.0, 5.0, 4.0, 8.0, 3.0, 10.0, 6.0, 3.0, 6.0, 7.0, 3.0, 7.0, 5.0, 7.0, 7.0, 6.0, 8.0, 4.0, 10.0, 2.0, 10.0, 5.0, 7.0, 3.0, 10.0, 7.0, 2.0, 3.0, 8.0, 2.0, 7.0, 6.0, 5.0, 9.0, 10.0, 3.0, 7.0, 3.0, 8.0, 9.0, 9.0, 10.0, 1.0, 4.0, 9.0, 6.0, 4.0, 5.0, 6.0, 9.0, 4.0, 1.0, 7.0, 8.0, 9.0, 9.0, 9.0, 3.0, 6.0, 2.0, 5.0, 2.0, 6.0, 1.0, 7.0, 5.0, 5.0, 5.0, 8.0, 6.0, 2.0, 9.0, 10.0, 9.0, 3.0, 8.0, 2.0, 5.0, 6.0, 8.0, 7.0, 9.0, 8.0, 6.0, 2.0, 7.0, 1.0, 3.0, 10.0, 1.0, 5.0, 9.0, 7.0, 8.0, 4.0, 2.0, 8.0, 3.0, 4.0, 1.0, 1.0, 4.0, 4.0, 4.0, 5.0, 1.0, 5.0, 2.0, 10.0, 3.0, 2.0, 1.0, 4.0, 7.0, 1.0, 9.0, 4.0, 4.0, 6.0, 4.0, 8.0, 6.0, 1.0, 6.0, 2.0, 6.0, 2.0, 1.0, 7.0, 6.0, 6.0, 2.0, 10.0]
global b_y = 10
global p = [0.076, 0.86, 0.815, 0.905, 0.146, 0.547, 0.322, 0.122, 0.359, 0.995, 0.562, 0.937, 0.062, 0.773, 0.657, 0.833, 0.518, 0.788, 0.05, 0.25, 0.589, 0.561, 0.322, 0.657, 0.248, 0.62, 0.246, 0.69, 0.164, 0.151, 0.733, 0.56, 0.155, 0.356, 0.186, 0.865, 0.769, 0.894, 0.49, 0.183, 0.405, 0.087, 0.062, 0.368, 0.792, 0.745, 0.93, 0.912, 0.854, 0.995, 0.397, 0.439, 0.234, 0.356, 0.713, 0.822, 0.056, 0.169, 0.418, 0.908, 0.75, 0.359, 0.578, 0.051, 0.546, 0.621, 0.429, 0.887, 0.723, 0.649, 0.275, 0.468, 0.166, 0.64, 0.323, 0.987, 0.992, 0.742, 0.87, 0.653, 0.606, 0.964, 0.795, 0.227, 0.514, 0.478, 0.255, 0.583, 0.494, 0.237, 0.905, 0.264, 0.076, 0.369, 0.225, 0.634, 0.269, 0.778, 0.553, 0.589, 0.527, 0.422, 0.904, 0.249, 0.094, 0.028, 0.037, 0.277, 0.952, 0.528, 0.993, 0.62, 0.357, 0.342, 0.457, 0.597, 0.534, 0.058, 0.775, 0.9, 0.827, 0.698, 0.351, 0.273, 0.674, 0.254, 0.784, 0.113, 0.4, 0.597, 0.343, 0.63, 0.895, 0.94, 0.3, 0.942, 0.175, 0.162, 0.45, 0.26, 0.122, 0.001, 0.198, 0.162, 0.831, 0.259, 0.353, 0.531, 0.253, 0.767, 0.307, 0.898, 0.128, 0.085, 0.779, 0.851, 0.512, 0.782, 0.022, 0.683, 0.892, 0.095, 0.991, 0.323, 0.345, 0.605, 0.148, 0.924, 0.287, 0.648, 0.359, 0.325, 0.369, 0.687, 0.166, 0.742, 0.265, 0.75, 0.889, 0.754, 0.664, 0.806, 0.062, 0.051, 0.524, 0.861, 0.119, 0.495, 0.927, 0.909, 0.467, 0.59, 0.02, 0.318, 0.893, 0.16, 0.5, 0.474, 0.989, 0.587, 0.49, 0.46, 0.704, 0.443, 0.709, 0.552, 0.858, 0.113, 0.123, 0.928, 0.055, 0.037, 0.827, 0.295, 0.544, 0.527, 0.656, 0.908, 0.186, 0.543, 0.024, 0.865, 0.079, 0.031, 0.262, 0.591, 0.876, 0.843, 0.442, 0.6, 0.894, 0.619, 0.256, 0.497, 0.947, 0.981, 0.212, 0.389]
global q = [0.967, 0.885, 0.817, 0.989, 0.947, 0.988, 0.512, 0.922, 0.413, 0.997, 0.795, 0.948, 0.564, 0.869, 0.887, 0.89, 0.991, 0.841, 0.935, 0.4, 0.89, 0.895, 0.416, 0.812, 0.674, 0.866, 0.47, 0.954, 0.558, 0.443, 0.749, 0.743, 0.593, 0.749, 0.407, 0.958, 0.809, 0.986, 0.737, 0.336, 0.87, 0.748, 0.817, 0.768, 0.822, 0.855, 0.996, 0.977, 0.864, 0.999, 0.806, 0.666, 0.464, 0.39, 0.878, 0.913, 0.59, 0.76, 0.58, 0.933, 0.818, 0.954, 0.807, 0.657, 0.55, 0.777, 0.859, 0.936, 0.876, 0.863, 0.947, 0.655, 0.305, 0.949, 0.727, 0.99, 0.999, 0.918, 0.948, 0.936, 0.978, 0.996, 0.99, 0.896, 0.918, 0.951, 0.449, 0.842, 0.811, 0.34, 0.997, 0.464, 0.839, 0.732, 0.954, 0.887, 0.971, 0.789, 0.665, 0.613, 0.62, 0.489, 0.905, 0.974, 0.766, 0.131, 0.041, 0.664, 0.966, 0.688, 0.997, 0.715, 0.57, 0.631, 0.922, 0.741, 0.919, 0.066, 0.871, 0.911, 0.978, 0.955, 0.467, 0.787, 0.913, 0.459, 0.961, 0.601, 0.94, 0.72, 0.378, 0.947, 0.944, 0.997, 0.617, 0.949, 0.291, 0.348, 0.685, 0.983, 0.37, 0.002, 0.646, 0.824, 0.874, 0.274, 0.475, 0.881, 0.369, 0.959, 0.986, 0.976, 0.33, 0.544, 0.826, 0.91, 0.652, 0.984, 0.844, 0.768, 0.973, 0.828, 0.998, 0.803, 0.654, 0.999, 0.497, 0.937, 0.837, 0.96, 0.567, 0.903, 0.456, 0.798, 0.23, 0.933, 0.906, 0.845, 0.997, 0.854, 0.953, 0.843, 0.843, 0.337, 0.772, 0.889, 0.814, 0.861, 0.977, 0.912, 0.65, 0.841, 0.084, 0.757, 0.977, 0.474, 0.559, 0.589, 0.989, 0.962, 0.533, 0.835, 0.814, 0.897, 0.74, 0.778, 0.916, 0.48, 0.711, 0.957, 0.861, 0.828, 0.83, 0.994, 0.944, 0.755, 0.94, 0.976, 0.455, 0.607, 0.912, 0.943, 0.236, 0.108, 0.338, 0.956, 0.905, 0.882, 0.78, 0.939, 0.978, 0.911, 0.724, 0.964, 0.999, 0.992, 0.785, 0.899]
global origin = 1
global destination = 50