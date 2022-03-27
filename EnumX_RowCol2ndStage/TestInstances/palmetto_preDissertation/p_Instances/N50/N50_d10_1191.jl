global arcs = [1 5; 1 34; 2 6; 2 19; 3 29; 4 5; 4 9; 4 17; 4 32; 4 36; 4 40; 5 19; 5 26; 5 29; 5 33; 5 38; 6 10; 6 11; 6 15; 6 18; 6 20; 6 32; 6 38; 6 44; 7 6; 7 12; 7 19; 7 20; 7 29; 7 33; 7 38; 8 9; 8 10; 8 13; 8 20; 8 22; 9 36; 9 47; 10 31; 10 34; 10 39; 10 47; 10 48; 11 5; 11 10; 11 27; 11 35; 11 38; 12 14; 12 25; 12 31; 12 34; 12 40; 12 41; 12 50; 13 26; 13 46; 13 48; 13 49; 14 2; 14 5; 14 10; 14 13; 14 18; 14 31; 14 38; 14 42; 14 48; 15 21; 15 30; 15 31; 15 32; 15 38; 15 40; 15 41; 15 42; 15 43; 16 13; 16 28; 16 33; 16 34; 16 43; 17 18; 17 40; 17 45; 18 20; 18 29; 18 32; 18 41; 19 3; 19 13; 19 17; 19 40; 20 16; 20 19; 20 21; 20 29; 20 31; 20 34; 20 35; 20 43; 20 49; 21 2; 21 6; 21 10; 21 30; 21 39; 22 4; 22 5; 22 10; 22 11; 22 14; 22 47; 22 50; 23 18; 23 27; 23 31; 23 34; 23 41; 24 8; 24 21; 24 34; 25 4; 25 20; 25 21; 25 22; 25 24; 25 38; 25 45; 26 11; 26 29; 26 33; 26 36; 26 37; 27 6; 27 13; 27 15; 27 17; 27 37; 27 39; 27 45; 27 50; 28 17; 28 31; 28 32; 29 6; 29 18; 29 31; 29 32; 29 48; 29 49; 30 7; 30 17; 30 22; 30 26; 30 33; 30 44; 31 5; 31 20; 31 41; 31 46; 31 49; 32 11; 32 18; 32 24; 32 28; 32 47; 33 15; 33 18; 33 20; 33 22; 33 35; 33 44; 34 10; 34 15; 34 16; 34 17; 34 30; 35 9; 35 13; 35 31; 35 41; 35 44; 36 11; 36 12; 36 15; 36 18; 36 28; 36 34; 36 44; 37 4; 37 7; 37 13; 37 15; 37 19; 37 33; 37 40; 37 43; 38 3; 38 17; 38 27; 38 41; 39 25; 39 28; 39 44; 39 46; 39 49; 40 6; 40 11; 40 19; 40 21; 40 24; 40 37; 40 44; 40 48; 41 2; 41 3; 41 33; 41 34; 41 36; 41 47; 42 11; 42 13; 42 49; 43 10; 43 14; 43 16; 43 18; 43 34; 44 9; 44 20; 44 26; 44 38; 45 2; 45 18; 45 41; 45 42; 46 23; 46 24; 46 31; 46 32; 47 2; 47 14; 47 16; 47 30; 47 46; 48 4; 48 11; 48 17; 48 19; 48 28; 48 36; 48 46; 49 3; 49 24; 49 30]
global d_x = [4.0, 3.0, 10.0, 8.0, 9.0, 1.0, 9.0, 3.0, 5.0, 7.0, 10.0, 6.0, 8.0, 6.0, 9.0, 8.0, 6.0, 8.0, 4.0, 7.0, 10.0, 2.0, 5.0, 3.0, 3.0, 9.0, 8.0, 6.0, 3.0, 6.0, 5.0, 3.0, 8.0, 6.0, 1.0, 4.0, 6.0, 5.0, 1.0, 8.0, 5.0, 5.0, 9.0, 3.0, 5.0, 10.0, 7.0, 6.0, 2.0, 7.0, 9.0, 5.0, 6.0, 6.0, 10.0, 7.0, 5.0, 6.0, 8.0, 7.0, 2.0, 10.0, 2.0, 1.0, 8.0, 2.0, 4.0, 9.0, 9.0, 9.0, 4.0, 1.0, 4.0, 8.0, 9.0, 7.0, 8.0, 9.0, 2.0, 1.0, 7.0, 3.0, 10.0, 6.0, 2.0, 6.0, 3.0, 5.0, 8.0, 6.0, 1.0, 6.0, 5.0, 7.0, 5.0, 2.0, 7.0, 5.0, 8.0, 3.0, 7.0, 3.0, 3.0, 10.0, 7.0, 2.0, 5.0, 5.0, 5.0, 3.0, 10.0, 4.0, 4.0, 6.0, 10.0, 8.0, 2.0, 6.0, 7.0, 9.0, 1.0, 2.0, 6.0, 2.0, 2.0, 6.0, 2.0, 3.0, 1.0, 8.0, 4.0, 10.0, 5.0, 3.0, 3.0, 5.0, 3.0, 3.0, 10.0, 8.0, 1.0, 2.0, 9.0, 9.0, 8.0, 7.0, 4.0, 5.0, 9.0, 9.0, 7.0, 4.0, 3.0, 10.0, 6.0, 4.0, 5.0, 6.0, 7.0, 6.0, 7.0, 9.0, 6.0, 6.0, 4.0, 4.0, 7.0, 5.0, 8.0, 8.0, 9.0, 5.0, 5.0, 2.0, 6.0, 5.0, 8.0, 3.0, 2.0, 2.0, 8.0, 10.0, 9.0, 4.0, 10.0, 1.0, 5.0, 9.0, 9.0, 7.0, 2.0, 3.0, 7.0, 4.0, 8.0, 2.0, 5.0, 1.0, 5.0, 2.0, 8.0, 5.0, 7.0, 7.0, 4.0, 10.0, 7.0, 10.0, 7.0, 8.0, 10.0, 8.0, 4.0, 3.0, 6.0, 3.0, 4.0, 4.0, 5.0, 4.0, 4.0, 5.0, 5.0, 1.0, 2.0, 9.0, 1.0, 3.0, 9.0, 2.0, 10.0, 8.0, 4.0, 10.0, 5.0, 1.0, 4.0, 3.0, 4.0, 6.0, 7.0, 7.0, 4.0, 6.0, 8.0, 10.0, 7.0, 8.0, 5.0, 4.0, 10.0, 3.0, 9.0, 6.0, 6.0, 4.0]
global b_x = 5
global d_y = [10.0, 6.0, 6.0, 9.0, 2.0, 8.0, 6.0, 8.0, 2.0, 2.0, 8.0, 10.0, 7.0, 1.0, 5.0, 7.0, 2.0, 9.0, 7.0, 1.0, 5.0, 1.0, 1.0, 5.0, 1.0, 9.0, 2.0, 4.0, 5.0, 7.0, 6.0, 8.0, 5.0, 10.0, 2.0, 4.0, 2.0, 6.0, 1.0, 3.0, 1.0, 9.0, 7.0, 6.0, 3.0, 8.0, 5.0, 3.0, 5.0, 8.0, 9.0, 2.0, 10.0, 2.0, 6.0, 2.0, 5.0, 8.0, 5.0, 7.0, 1.0, 4.0, 6.0, 9.0, 8.0, 10.0, 2.0, 3.0, 2.0, 8.0, 8.0, 7.0, 9.0, 2.0, 7.0, 7.0, 7.0, 7.0, 9.0, 6.0, 1.0, 2.0, 7.0, 7.0, 10.0, 9.0, 2.0, 7.0, 2.0, 9.0, 5.0, 1.0, 5.0, 7.0, 1.0, 10.0, 9.0, 4.0, 6.0, 9.0, 7.0, 2.0, 10.0, 8.0, 7.0, 8.0, 3.0, 5.0, 10.0, 1.0, 8.0, 2.0, 6.0, 6.0, 6.0, 4.0, 9.0, 10.0, 3.0, 2.0, 10.0, 5.0, 4.0, 8.0, 3.0, 7.0, 4.0, 10.0, 9.0, 6.0, 8.0, 4.0, 2.0, 5.0, 6.0, 9.0, 9.0, 8.0, 7.0, 5.0, 4.0, 10.0, 7.0, 9.0, 3.0, 3.0, 6.0, 6.0, 3.0, 4.0, 7.0, 5.0, 7.0, 10.0, 2.0, 2.0, 5.0, 2.0, 7.0, 3.0, 6.0, 4.0, 10.0, 6.0, 7.0, 5.0, 1.0, 9.0, 9.0, 3.0, 5.0, 8.0, 6.0, 9.0, 1.0, 7.0, 3.0, 9.0, 10.0, 10.0, 1.0, 2.0, 5.0, 10.0, 1.0, 8.0, 3.0, 7.0, 2.0, 1.0, 5.0, 9.0, 3.0, 4.0, 2.0, 3.0, 4.0, 4.0, 6.0, 4.0, 4.0, 8.0, 6.0, 3.0, 5.0, 5.0, 8.0, 5.0, 2.0, 7.0, 7.0, 2.0, 5.0, 4.0, 9.0, 3.0, 4.0, 1.0, 5.0, 7.0, 5.0, 6.0, 4.0, 4.0, 4.0, 4.0, 1.0, 10.0, 5.0, 4.0, 8.0, 8.0, 9.0, 2.0, 5.0, 3.0, 5.0, 9.0, 6.0, 2.0, 6.0, 6.0, 3.0, 3.0, 3.0, 6.0, 1.0, 10.0, 8.0, 2.0, 4.0, 3.0, 6.0, 9.0, 3.0, 10.0]
global b_y = 10
global p = [0.832, 0.704, 0.623, 0.056, 0.961, 0.293, 0.429, 0.414, 0.089, 0.29, 0.247, 0.271, 0.78, 0.296, 0.931, 0.009, 0.214, 0.919, 0.433, 0.257, 0.456, 0.867, 0.605, 0.179, 0.195, 0.057, 0.047, 0.282, 0.335, 0.353, 0.187, 0.923, 0.228, 0.532, 0.04, 0.698, 0.193, 0.766, 0.544, 0.513, 0.544, 0.518, 0.54, 0.186, 0.534, 0.393, 0.568, 0.251, 0.847, 0.956, 0.08, 0.455, 0.985, 0.226, 0.117, 0.748, 0.367, 0.555, 0.886, 0.016, 0.751, 0.278, 0.344, 0.285, 0.276, 0.321, 0.851, 0.623, 0.844, 0.829, 0.893, 0.632, 0.494, 0.344, 0.069, 0.419, 0.203, 0.775, 0.281, 0.689, 0.408, 0.822, 0.464, 0.623, 0.185, 0.793, 0.491, 0.424, 0.751, 0.275, 0.563, 0.405, 0.444, 0.726, 0.61, 0.968, 0.923, 0.64, 0.971, 0.718, 0.343, 0.208, 0.81, 0.833, 0.791, 0.345, 0.512, 0.035, 0.232, 0.08, 0.785, 0.12, 0.36, 0.456, 0.923, 0.876, 0.976, 0.756, 0.969, 0.631, 0.057, 0.792, 0.889, 0.957, 0.125, 0.644, 0.498, 0.803, 0.97, 0.964, 0.868, 0.364, 0.014, 0.804, 0.468, 0.46, 0.487, 0.745, 0.675, 0.84, 0.07, 0.952, 0.66, 0.037, 0.521, 0.956, 0.697, 0.367, 0.17, 0.046, 0.601, 0.926, 0.605, 0.234, 0.992, 0.962, 0.797, 0.653, 0.285, 0.379, 0.612, 0.55, 0.089, 0.085, 0.601, 0.091, 0.566, 0.969, 0.621, 0.559, 0.91, 0.677, 0.937, 0.947, 0.172, 0.789, 0.024, 0.306, 0.956, 0.222, 0.386, 0.645, 0.403, 0.667, 0.397, 0.678, 0.74, 0.81, 0.495, 0.05, 0.033, 0.575, 0.413, 0.804, 0.474, 0.616, 0.728, 0.506, 0.831, 0.827, 0.386, 0.024, 0.266, 0.64, 0.873, 0.135, 0.892, 0.079, 0.925, 0.845, 0.286, 0.687, 0.944, 0.47, 0.193, 0.493, 0.265, 0.453, 0.699, 0.034, 0.118, 0.327, 0.991, 0.234, 0.9, 0.661, 0.397, 0.818, 0.561, 0.596, 0.766, 0.622, 0.873, 0.03, 0.743, 0.405, 0.682, 0.911, 0.023, 0.71, 0.888, 0.213, 0.61, 0.655, 0.009, 0.021, 0.294, 0.285, 0.25, 0.056, 0.189, 0.314, 0.669, 0.2, 0.429, 0.459]
global q = [0.929, 0.753, 0.795, 0.253, 0.972, 0.652, 0.549, 0.696, 0.878, 0.967, 0.725, 0.377, 0.928, 0.755, 0.948, 0.567, 0.631, 0.956, 0.768, 0.257, 0.584, 0.984, 0.776, 0.204, 0.441, 0.909, 0.249, 0.976, 0.795, 0.568, 0.707, 0.994, 0.554, 0.541, 0.629, 0.985, 0.531, 0.914, 0.547, 0.841, 0.707, 0.858, 0.616, 0.573, 0.791, 0.752, 0.82, 0.535, 0.871, 0.991, 0.6, 0.982, 0.985, 0.88, 0.636, 0.886, 0.933, 0.834, 0.967, 0.664, 0.891, 0.547, 0.44, 0.842, 0.496, 0.412, 0.91, 0.728, 0.907, 0.853, 0.898, 0.887, 0.605, 0.524, 0.897, 0.728, 0.78, 0.8, 0.934, 0.82, 0.915, 0.989, 0.475, 0.709, 0.212, 0.915, 0.949, 0.605, 0.858, 0.311, 0.601, 0.511, 0.484, 0.918, 0.847, 0.985, 0.957, 0.878, 0.974, 0.892, 0.427, 0.796, 0.962, 0.965, 0.93, 0.519, 0.845, 0.289, 0.714, 0.198, 0.834, 0.501, 0.584, 0.656, 0.962, 0.99, 0.98, 0.91, 0.989, 0.642, 0.438, 0.93, 0.89, 0.999, 0.918, 0.848, 0.505, 0.97, 0.978, 0.992, 0.933, 0.521, 0.132, 0.964, 0.782, 0.667, 0.534, 0.757, 0.815, 0.862, 0.569, 0.981, 0.724, 0.804, 0.583, 0.978, 0.8, 0.531, 0.967, 0.649, 0.81, 0.944, 0.732, 0.786, 0.999, 0.983, 0.842, 0.915, 0.757, 0.72, 0.832, 0.701, 0.829, 0.393, 0.676, 0.626, 0.764, 0.998, 0.778, 0.61, 0.992, 0.701, 0.989, 0.972, 0.69, 0.971, 0.152, 0.923, 0.956, 0.468, 0.664, 0.966, 0.489, 0.769, 0.517, 0.696, 0.874, 0.997, 0.768, 0.987, 0.343, 0.828, 0.599, 0.927, 0.536, 0.663, 0.948, 0.989, 0.95, 0.996, 0.643, 0.64, 0.283, 0.776, 0.918, 0.924, 0.98, 0.582, 0.935, 0.847, 0.758, 0.797, 0.962, 0.649, 0.725, 0.522, 0.892, 0.573, 0.943, 0.501, 0.352, 0.416, 0.997, 0.826, 0.922, 0.7, 0.833, 0.825, 0.779, 0.974, 0.909, 0.91, 0.881, 0.921, 0.784, 0.484, 0.696, 0.988, 0.331, 0.953, 0.895, 0.82, 0.793, 0.809, 0.158, 0.165, 0.964, 0.626, 0.578, 0.977, 0.569, 0.87, 0.675, 0.879, 0.898, 0.978]
global origin = 1
global destination = 50