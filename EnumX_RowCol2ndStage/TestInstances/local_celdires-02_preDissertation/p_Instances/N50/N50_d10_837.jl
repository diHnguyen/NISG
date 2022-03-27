global arcs = [1 44; 2 13; 2 23; 2 25; 3 28; 3 33; 3 39; 3 42; 3 45; 3 49; 4 34; 4 41; 4 44; 5 36; 5 42; 5 48; 6 5; 6 25; 6 30; 6 31; 6 32; 6 39; 6 47; 6 48; 6 49; 7 15; 7 30; 7 41; 7 45; 8 12; 8 23; 8 38; 8 47; 9 6; 9 7; 9 16; 9 19; 9 25; 9 33; 9 35; 9 38; 9 41; 10 11; 10 24; 10 38; 10 50; 11 25; 11 29; 11 37; 11 38; 12 13; 12 17; 12 33; 12 45; 13 33; 13 39; 13 40; 14 4; 14 8; 14 21; 14 22; 14 25; 14 34; 14 40; 15 24; 15 38; 15 47; 16 7; 16 18; 17 4; 17 5; 17 7; 17 10; 17 22; 17 34; 17 38; 17 48; 18 6; 18 14; 18 20; 18 28; 19 3; 19 6; 19 8; 19 20; 19 28; 19 37; 20 17; 20 39; 20 50; 21 3; 21 39; 21 43; 21 46; 22 5; 22 6; 22 9; 22 11; 22 13; 22 21; 22 37; 22 39; 22 41; 22 44; 23 4; 23 7; 23 13; 23 17; 23 24; 23 25; 23 43; 23 47; 24 26; 24 36; 24 42; 24 43; 24 45; 25 5; 25 13; 25 22; 25 43; 26 6; 26 8; 26 12; 26 21; 26 47; 27 20; 27 23; 27 25; 27 36; 27 45; 28 8; 28 14; 28 24; 28 35; 28 40; 28 41; 28 48; 29 2; 29 16; 29 30; 29 32; 29 35; 29 38; 29 42; 29 43; 29 45; 30 20; 30 36; 30 37; 30 40; 30 48; 30 49; 31 16; 31 19; 31 27; 32 4; 32 7; 32 8; 32 12; 32 36; 32 39; 32 41; 33 7; 33 22; 33 32; 33 34; 33 38; 33 43; 33 47; 33 48; 34 2; 34 6; 34 7; 34 26; 34 37; 34 38; 35 9; 35 44; 35 49; 36 5; 36 8; 36 11; 36 29; 36 47; 37 3; 37 8; 37 9; 37 16; 37 17; 37 34; 37 46; 37 47; 38 6; 38 17; 38 19; 38 20; 38 21; 38 26; 38 46; 38 48; 39 6; 39 14; 39 20; 39 32; 39 35; 39 38; 40 15; 40 19; 40 24; 40 25; 40 29; 40 33; 40 39; 41 6; 41 16; 41 29; 41 45; 42 7; 42 10; 42 15; 42 36; 42 48; 42 50; 43 10; 43 11; 43 13; 43 21; 43 38; 43 50; 44 5; 44 6; 44 26; 44 29; 44 32; 44 39; 45 4; 45 8; 45 10; 45 12; 45 28; 45 33; 45 44; 46 2; 46 6; 46 16; 46 24; 46 32; 46 45; 46 47; 47 10; 47 11; 47 27; 47 29; 48 20; 48 21; 48 27; 48 28; 48 30; 48 40; 49 4; 49 14; 49 22; 49 32; 49 40]
global d_x = [10.0, 3.0, 5.0, 1.0, 1.0, 1.0, 1.0, 4.0, 10.0, 1.0, 7.0, 4.0, 2.0, 1.0, 4.0, 5.0, 7.0, 7.0, 8.0, 3.0, 2.0, 10.0, 4.0, 3.0, 4.0, 4.0, 5.0, 3.0, 7.0, 8.0, 9.0, 1.0, 7.0, 3.0, 3.0, 9.0, 1.0, 9.0, 3.0, 10.0, 1.0, 2.0, 7.0, 9.0, 4.0, 6.0, 8.0, 9.0, 8.0, 5.0, 6.0, 8.0, 9.0, 1.0, 10.0, 1.0, 8.0, 3.0, 2.0, 5.0, 7.0, 10.0, 2.0, 2.0, 4.0, 7.0, 4.0, 1.0, 6.0, 2.0, 10.0, 4.0, 10.0, 4.0, 9.0, 8.0, 7.0, 5.0, 8.0, 8.0, 6.0, 2.0, 6.0, 10.0, 7.0, 7.0, 3.0, 4.0, 3.0, 1.0, 6.0, 6.0, 9.0, 7.0, 3.0, 1.0, 7.0, 8.0, 5.0, 3.0, 7.0, 2.0, 2.0, 1.0, 1.0, 8.0, 10.0, 3.0, 8.0, 10.0, 4.0, 2.0, 5.0, 8.0, 6.0, 1.0, 7.0, 2.0, 7.0, 6.0, 6.0, 8.0, 5.0, 9.0, 10.0, 1.0, 5.0, 5.0, 9.0, 4.0, 3.0, 3.0, 2.0, 8.0, 3.0, 9.0, 10.0, 10.0, 9.0, 2.0, 6.0, 6.0, 6.0, 1.0, 1.0, 5.0, 6.0, 1.0, 3.0, 5.0, 6.0, 1.0, 1.0, 3.0, 3.0, 2.0, 5.0, 10.0, 2.0, 4.0, 1.0, 6.0, 10.0, 6.0, 6.0, 10.0, 3.0, 2.0, 10.0, 4.0, 8.0, 8.0, 10.0, 10.0, 5.0, 4.0, 10.0, 2.0, 4.0, 7.0, 10.0, 9.0, 4.0, 1.0, 6.0, 7.0, 7.0, 3.0, 10.0, 7.0, 7.0, 3.0, 3.0, 10.0, 6.0, 4.0, 2.0, 8.0, 6.0, 7.0, 3.0, 8.0, 8.0, 7.0, 8.0, 2.0, 9.0, 10.0, 3.0, 4.0, 9.0, 9.0, 1.0, 9.0, 7.0, 10.0, 3.0, 10.0, 4.0, 5.0, 9.0, 6.0, 9.0, 2.0, 6.0, 2.0, 5.0, 5.0, 5.0, 4.0, 4.0, 5.0, 3.0, 9.0, 10.0, 3.0, 10.0, 8.0, 6.0, 3.0, 2.0, 10.0, 5.0, 1.0, 2.0, 2.0, 3.0, 8.0, 7.0, 9.0, 9.0, 4.0, 7.0, 7.0, 4.0, 9.0, 10.0, 10.0, 8.0, 2.0, 3.0, 9.0, 9.0, 1.0, 8.0]
global b_x = 5
global d_y = [4.0, 8.0, 6.0, 9.0, 5.0, 9.0, 2.0, 3.0, 8.0, 1.0, 1.0, 4.0, 10.0, 10.0, 1.0, 4.0, 1.0, 1.0, 8.0, 4.0, 4.0, 7.0, 3.0, 2.0, 6.0, 7.0, 1.0, 4.0, 2.0, 4.0, 2.0, 10.0, 8.0, 9.0, 6.0, 2.0, 4.0, 4.0, 7.0, 7.0, 9.0, 2.0, 4.0, 7.0, 6.0, 2.0, 7.0, 2.0, 10.0, 9.0, 3.0, 5.0, 5.0, 4.0, 2.0, 7.0, 9.0, 10.0, 4.0, 2.0, 10.0, 1.0, 3.0, 7.0, 5.0, 5.0, 5.0, 5.0, 7.0, 6.0, 10.0, 9.0, 1.0, 2.0, 1.0, 2.0, 10.0, 5.0, 2.0, 5.0, 8.0, 4.0, 2.0, 5.0, 6.0, 6.0, 1.0, 3.0, 2.0, 5.0, 6.0, 3.0, 9.0, 6.0, 7.0, 3.0, 2.0, 1.0, 6.0, 9.0, 10.0, 9.0, 9.0, 8.0, 9.0, 4.0, 2.0, 2.0, 4.0, 5.0, 3.0, 4.0, 7.0, 2.0, 5.0, 8.0, 2.0, 4.0, 4.0, 4.0, 8.0, 5.0, 1.0, 8.0, 4.0, 8.0, 6.0, 4.0, 5.0, 2.0, 7.0, 7.0, 4.0, 9.0, 6.0, 2.0, 9.0, 4.0, 7.0, 2.0, 3.0, 7.0, 6.0, 6.0, 8.0, 1.0, 6.0, 9.0, 8.0, 7.0, 6.0, 8.0, 5.0, 5.0, 3.0, 1.0, 3.0, 10.0, 10.0, 3.0, 10.0, 5.0, 5.0, 6.0, 2.0, 2.0, 4.0, 2.0, 5.0, 4.0, 2.0, 8.0, 8.0, 1.0, 1.0, 7.0, 2.0, 5.0, 3.0, 4.0, 9.0, 3.0, 4.0, 3.0, 5.0, 2.0, 9.0, 6.0, 7.0, 9.0, 5.0, 7.0, 7.0, 10.0, 1.0, 8.0, 6.0, 3.0, 5.0, 3.0, 4.0, 7.0, 9.0, 7.0, 8.0, 1.0, 10.0, 6.0, 1.0, 6.0, 10.0, 7.0, 6.0, 10.0, 9.0, 10.0, 7.0, 6.0, 2.0, 6.0, 1.0, 2.0, 9.0, 3.0, 4.0, 10.0, 2.0, 2.0, 5.0, 10.0, 6.0, 2.0, 1.0, 1.0, 10.0, 10.0, 4.0, 3.0, 8.0, 7.0, 4.0, 3.0, 2.0, 9.0, 4.0, 3.0, 7.0, 9.0, 5.0, 7.0, 9.0, 5.0, 7.0, 10.0, 5.0, 10.0, 4.0, 8.0, 4.0, 5.0, 8.0, 3.0, 10.0, 5.0, 9.0]
global b_y = 10
global p = [0.844, 0.467, 0.208, 0.896, 0.462, 0.192, 0.717, 0.29, 0.924, 0.642, 0.406, 0.944, 0.529, 0.86, 0.39, 0.918, 0.427, 0.379, 0.829, 0.905, 0.629, 0.993, 0.737, 0.205, 0.001, 0.431, 0.582, 0.848, 0.874, 0.881, 0.988, 0.336, 0.098, 0.949, 0.629, 0.08, 0.836, 0.427, 0.868, 0.024, 0.531, 0.745, 0.237, 0.696, 0.509, 0.611, 0.213, 0.899, 0.879, 0.019, 0.805, 0.125, 0.138, 0.166, 0.369, 0.408, 0.851, 0.038, 0.166, 0.734, 0.18, 0.593, 0.533, 0.472, 0.862, 0.843, 0.091, 0.397, 0.422, 0.793, 0.983, 0.307, 0.229, 0.831, 0.308, 0.727, 0.655, 0.22, 0.416, 0.151, 0.551, 0.85, 0.953, 0.499, 0.782, 0.856, 0.138, 0.828, 0.962, 0.114, 0.49, 0.955, 0.424, 0.829, 0.149, 0.073, 0.478, 0.069, 0.745, 0.339, 0.488, 0.408, 0.977, 0.484, 0.639, 0.118, 0.046, 0.259, 0.785, 0.958, 0.416, 0.551, 0.455, 0.454, 0.579, 0.574, 0.622, 0.363, 0.88, 0.229, 0.641, 0.071, 0.265, 0.151, 0.764, 0.867, 0.298, 0.565, 0.504, 0.127, 0.091, 0.357, 0.832, 0.677, 0.196, 0.244, 0.207, 0.821, 0.327, 0.764, 0.266, 0.589, 0.037, 0.228, 0.328, 0.84, 0.985, 0.358, 0.459, 0.726, 0.721, 0.831, 0.395, 0.423, 0.752, 0.926, 0.276, 0.331, 0.333, 0.677, 0.773, 0.899, 0.899, 0.604, 0.388, 0.376, 0.994, 0.032, 0.785, 0.228, 0.087, 0.907, 0.682, 0.504, 0.963, 0.617, 0.042, 0.784, 0.613, 0.176, 0.5, 0.812, 0.642, 0.069, 0.57, 0.117, 0.541, 0.491, 0.289, 0.567, 0.461, 0.769, 0.796, 0.514, 0.259, 0.956, 0.328, 0.218, 0.381, 0.986, 0.118, 0.907, 0.308, 0.171, 0.444, 0.892, 0.89, 0.232, 0.58, 0.77, 0.744, 0.412, 0.309, 0.059, 0.845, 0.929, 0.123, 0.358, 0.88, 0.481, 0.011, 0.385, 0.895, 0.617, 0.981, 0.992, 0.311, 0.57, 0.719, 0.847, 0.465, 0.886, 0.037, 0.376, 0.783, 0.347, 0.852, 0.929, 0.804, 0.888, 0.89, 0.03, 0.114, 0.796, 0.348, 0.811, 0.691, 0.208, 0.069, 0.548, 0.729, 0.857, 0.321, 0.389, 0.431, 0.009, 0.625, 0.885, 0.829, 0.153, 0.193, 0.402, 0.537, 0.378, 0.211]
global q = [0.979, 0.701, 0.377, 0.932, 0.996, 0.346, 0.887, 0.945, 0.944, 0.727, 0.558, 0.998, 0.596, 0.967, 0.682, 0.935, 0.863, 0.886, 0.891, 0.914, 0.706, 0.995, 0.745, 0.456, 0.658, 0.438, 0.795, 0.877, 0.931, 0.936, 0.993, 0.656, 0.907, 0.977, 0.982, 0.55, 0.862, 0.744, 0.882, 0.281, 0.908, 0.89, 0.631, 0.78, 0.649, 0.781, 0.881, 0.981, 0.895, 0.678, 0.893, 0.362, 0.813, 0.711, 0.447, 0.471, 0.941, 0.573, 0.991, 0.854, 0.531, 0.691, 0.575, 0.628, 0.986, 0.976, 0.124, 0.635, 0.997, 0.904, 0.984, 0.723, 0.566, 0.883, 0.76, 0.803, 0.734, 0.316, 0.637, 0.556, 0.921, 0.928, 0.989, 0.769, 0.837, 0.934, 0.802, 0.88, 0.993, 0.636, 0.52, 0.956, 0.434, 0.843, 0.837, 0.431, 0.709, 0.263, 0.887, 0.75, 0.801, 0.535, 0.99, 0.917, 0.879, 0.19, 0.554, 0.799, 0.956, 0.998, 0.469, 0.585, 0.494, 0.684, 0.855, 0.828, 0.963, 0.756, 0.911, 0.437, 0.836, 0.605, 0.492, 0.309, 0.8, 0.92, 0.528, 0.919, 0.522, 0.989, 0.951, 0.507, 0.994, 0.918, 0.428, 0.554, 0.762, 0.946, 0.756, 0.949, 0.966, 0.693, 0.344, 0.576, 0.948, 0.846, 0.991, 0.729, 0.761, 0.897, 0.736, 0.875, 0.938, 0.455, 0.961, 0.996, 0.679, 0.542, 0.563, 0.872, 0.847, 0.943, 0.925, 0.954, 0.83, 0.69, 0.994, 0.309, 0.884, 0.989, 0.518, 0.925, 0.985, 0.943, 0.996, 0.841, 0.4, 0.99, 0.729, 0.96, 0.774, 0.914, 0.842, 0.816, 0.901, 0.285, 0.797, 0.717, 0.951, 0.732, 0.66, 0.981, 0.825, 0.767, 0.354, 0.963, 0.727, 0.223, 0.93, 0.996, 0.139, 0.959, 0.914, 0.767, 0.79, 0.957, 0.911, 0.674, 0.678, 0.983, 0.849, 0.884, 0.558, 0.076, 0.855, 0.952, 0.604, 0.409, 0.905, 0.972, 0.952, 0.885, 0.938, 0.818, 0.996, 0.998, 0.685, 0.897, 0.774, 0.987, 0.838, 0.931, 0.699, 0.867, 0.999, 0.486, 0.854, 0.938, 0.88, 0.994, 0.925, 0.606, 0.536, 0.963, 0.963, 0.923, 0.792, 0.431, 0.67, 0.756, 0.776, 0.977, 0.984, 0.455, 0.966, 0.396, 0.677, 0.949, 0.887, 0.167, 0.544, 0.609, 0.925, 0.623, 0.282]
global origin = 1
global destination = 50