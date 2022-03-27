global arcs = [1 3; 1 9; 1 14; 1 19; 1 41; 1 47; 2 34; 2 43; 2 44; 3 9; 3 22; 3 38; 3 42; 3 44; 4 15; 4 16; 5 10; 5 11; 5 15; 5 23; 5 37; 5 41; 6 18; 6 22; 6 25; 6 40; 7 5; 7 6; 7 9; 7 15; 7 24; 7 27; 7 49; 8 27; 8 40; 9 5; 9 19; 9 21; 9 31; 9 37; 10 12; 10 27; 11 7; 11 12; 11 31; 12 29; 12 49; 13 21; 13 23; 13 31; 14 13; 14 26; 14 35; 14 36; 14 45; 15 3; 15 29; 15 47; 16 15; 16 36; 16 44; 16 46; 17 9; 17 19; 17 27; 17 35; 17 43; 17 44; 18 30; 18 32; 19 14; 19 16; 19 27; 19 29; 19 38; 19 44; 20 10; 20 18; 20 22; 20 33; 21 4; 21 18; 21 28; 21 34; 21 47; 21 49; 22 3; 22 25; 22 32; 22 34; 22 45; 22 47; 22 50; 23 4; 23 17; 23 25; 24 7; 24 10; 24 16; 25 8; 25 16; 25 20; 25 33; 25 47; 25 48; 26 9; 26 20; 26 22; 26 42; 26 45; 26 46; 27 8; 27 14; 27 16; 27 18; 27 28; 27 31; 28 9; 28 10; 28 21; 28 34; 28 39; 28 42; 29 10; 29 19; 29 22; 29 32; 29 40; 30 12; 30 15; 30 20; 30 37; 30 48; 31 4; 31 12; 31 16; 31 19; 31 34; 31 35; 31 39; 32 4; 32 5; 32 6; 32 7; 32 10; 32 15; 32 17; 32 35; 32 41; 33 24; 33 35; 33 44; 33 47; 33 50; 34 6; 34 11; 34 30; 34 33; 34 49; 35 11; 35 18; 35 30; 35 31; 35 42; 36 23; 36 25; 37 16; 37 48; 38 4; 38 6; 38 14; 38 17; 38 19; 38 36; 39 21; 39 26; 39 46; 40 10; 40 18; 40 19; 40 22; 41 11; 41 18; 41 24; 41 35; 42 13; 43 18; 43 22; 43 25; 44 19; 44 23; 44 35; 45 2; 45 9; 45 14; 45 20; 45 33; 45 38; 45 39; 45 41; 46 8; 46 10; 46 15; 46 22; 46 27; 46 36; 46 43; 47 4; 47 10; 47 15; 47 44; 48 4; 48 17; 48 29; 48 37; 48 42; 48 45; 49 11; 49 36; 49 37; 49 42; 49 44; 49 46; 49 47; 49 48]
global d_x = [3.0, 4.0, 3.0, 5.0, 8.0, 8.0, 10.0, 6.0, 2.0, 1.0, 1.0, 9.0, 10.0, 5.0, 8.0, 6.0, 5.0, 4.0, 9.0, 4.0, 7.0, 5.0, 9.0, 4.0, 4.0, 1.0, 10.0, 5.0, 3.0, 10.0, 5.0, 5.0, 5.0, 5.0, 6.0, 10.0, 1.0, 4.0, 2.0, 3.0, 10.0, 5.0, 8.0, 2.0, 2.0, 6.0, 10.0, 5.0, 9.0, 6.0, 10.0, 5.0, 4.0, 3.0, 9.0, 8.0, 10.0, 1.0, 8.0, 5.0, 9.0, 8.0, 10.0, 4.0, 7.0, 2.0, 7.0, 10.0, 5.0, 9.0, 2.0, 6.0, 3.0, 5.0, 5.0, 3.0, 10.0, 7.0, 1.0, 4.0, 2.0, 1.0, 1.0, 1.0, 4.0, 8.0, 2.0, 9.0, 8.0, 4.0, 2.0, 2.0, 2.0, 2.0, 8.0, 8.0, 1.0, 8.0, 1.0, 10.0, 2.0, 2.0, 4.0, 5.0, 8.0, 2.0, 1.0, 3.0, 5.0, 3.0, 8.0, 6.0, 9.0, 3.0, 8.0, 3.0, 1.0, 10.0, 9.0, 8.0, 4.0, 6.0, 7.0, 5.0, 8.0, 6.0, 3.0, 9.0, 9.0, 10.0, 10.0, 9.0, 8.0, 10.0, 2.0, 6.0, 8.0, 6.0, 2.0, 7.0, 10.0, 4.0, 9.0, 8.0, 9.0, 2.0, 4.0, 7.0, 1.0, 8.0, 10.0, 1.0, 9.0, 5.0, 10.0, 2.0, 1.0, 10.0, 5.0, 1.0, 6.0, 4.0, 8.0, 3.0, 5.0, 4.0, 4.0, 3.0, 7.0, 5.0, 6.0, 10.0, 4.0, 4.0, 1.0, 10.0, 5.0, 5.0, 8.0, 5.0, 6.0, 9.0, 1.0, 1.0, 10.0, 1.0, 10.0, 8.0, 3.0, 3.0, 9.0, 8.0, 6.0, 4.0, 7.0, 2.0, 10.0, 5.0, 2.0, 2.0, 9.0, 3.0, 7.0, 6.0, 7.0, 6.0, 2.0, 9.0, 4.0, 1.0, 9.0, 6.0, 9.0, 6.0, 9.0, 9.0, 10.0, 6.0, 1.0, 4.0, 2.0, 2.0, 7.0, 9.0, 5.0]
global b_x = 5
global d_y = [8.0, 1.0, 2.0, 2.0, 10.0, 2.0, 6.0, 2.0, 3.0, 4.0, 6.0, 5.0, 4.0, 7.0, 5.0, 2.0, 2.0, 10.0, 7.0, 5.0, 3.0, 10.0, 10.0, 10.0, 7.0, 5.0, 1.0, 6.0, 2.0, 9.0, 10.0, 5.0, 8.0, 9.0, 1.0, 7.0, 10.0, 5.0, 1.0, 3.0, 9.0, 5.0, 3.0, 8.0, 2.0, 8.0, 9.0, 10.0, 7.0, 7.0, 6.0, 7.0, 9.0, 7.0, 7.0, 2.0, 8.0, 6.0, 3.0, 4.0, 10.0, 5.0, 5.0, 1.0, 5.0, 7.0, 5.0, 5.0, 3.0, 1.0, 2.0, 8.0, 6.0, 5.0, 2.0, 5.0, 3.0, 3.0, 5.0, 10.0, 1.0, 10.0, 7.0, 4.0, 3.0, 8.0, 1.0, 5.0, 8.0, 3.0, 2.0, 9.0, 10.0, 6.0, 1.0, 3.0, 1.0, 9.0, 4.0, 10.0, 1.0, 4.0, 8.0, 10.0, 7.0, 4.0, 8.0, 7.0, 6.0, 6.0, 4.0, 5.0, 3.0, 5.0, 10.0, 5.0, 7.0, 10.0, 5.0, 5.0, 1.0, 1.0, 8.0, 9.0, 4.0, 2.0, 5.0, 3.0, 5.0, 9.0, 5.0, 4.0, 10.0, 2.0, 6.0, 10.0, 10.0, 3.0, 7.0, 2.0, 10.0, 8.0, 7.0, 5.0, 3.0, 9.0, 5.0, 2.0, 7.0, 4.0, 8.0, 8.0, 8.0, 7.0, 7.0, 4.0, 6.0, 10.0, 8.0, 7.0, 9.0, 10.0, 5.0, 6.0, 2.0, 2.0, 3.0, 2.0, 2.0, 8.0, 9.0, 1.0, 5.0, 8.0, 10.0, 5.0, 9.0, 9.0, 6.0, 4.0, 2.0, 2.0, 4.0, 4.0, 9.0, 4.0, 1.0, 3.0, 7.0, 3.0, 7.0, 6.0, 4.0, 1.0, 9.0, 10.0, 3.0, 2.0, 3.0, 2.0, 3.0, 4.0, 10.0, 1.0, 9.0, 3.0, 2.0, 4.0, 6.0, 8.0, 10.0, 5.0, 9.0, 4.0, 6.0, 4.0, 4.0, 4.0, 2.0, 3.0, 1.0, 9.0, 1.0, 6.0, 8.0]
global b_y = 10
global p = [0.668, 0.886, 0.513, 0.135, 0.516, 0.515, 0.387, 0.508, 0.778, 0.622, 0.511, 0.115, 0.379, 0.814, 0.076, 0.304, 0.867, 0.539, 0.937, 0.985, 0.826, 0.205, 0.475, 0.563, 0.546, 0.573, 0.61, 0.594, 0.13, 0.432, 0.792, 0.476, 0.716, 0.721, 0.245, 0.038, 0.878, 0.971, 0.277, 0.682, 0.59, 0.718, 0.332, 0.63, 0.153, 0.283, 0.53, 0.775, 0.982, 0.446, 0.594, 0.52, 0.801, 0.245, 0.951, 0.488, 0.231, 0.724, 0.645, 0.925, 0.949, 0.761, 0.964, 0.745, 0.659, 0.496, 0.262, 0.172, 0.768, 0.502, 0.959, 0.762, 0.59, 0.406, 0.304, 0.721, 0.593, 0.877, 0.875, 0.723, 0.046, 0.242, 0.18, 0.205, 0.742, 0.772, 0.139, 0.977, 0.825, 0.298, 0.451, 0.756, 0.516, 0.372, 0.929, 0.781, 0.543, 0.449, 0.215, 0.632, 0.788, 0.64, 0.245, 0.172, 0.429, 0.487, 0.067, 0.542, 0.521, 0.935, 0.708, 0.087, 0.034, 0.807, 0.092, 0.895, 0.824, 0.311, 0.722, 0.043, 0.936, 0.715, 0.698, 0.373, 0.126, 0.006, 0.759, 0.971, 0.284, 0.177, 0.032, 0.965, 0.532, 0.457, 0.309, 0.956, 0.885, 0.727, 0.019, 0.303, 0.864, 0.386, 0.435, 0.733, 0.608, 0.724, 0.212, 0.13, 0.1, 0.807, 0.043, 0.181, 0.75, 0.741, 0.588, 0.621, 0.55, 0.987, 0.675, 0.911, 0.398, 0.725, 0.761, 0.199, 0.24, 0.927, 0.114, 0.235, 0.712, 0.57, 0.218, 0.619, 0.646, 0.833, 0.458, 0.099, 0.664, 0.475, 0.852, 0.603, 0.547, 0.29, 0.142, 0.774, 0.897, 0.754, 0.06, 0.631, 0.23, 0.662, 0.472, 0.071, 0.909, 0.569, 0.618, 0.745, 0.494, 0.302, 0.438, 0.983, 0.588, 0.862, 0.583, 0.416, 0.508, 0.678, 0.75, 0.596, 0.774, 0.327, 0.339, 0.245, 0.614, 0.245, 0.999, 0.338, 0.714, 0.645, 0.723, 0.524, 0.12, 0.225, 0.838, 0.239, 0.517]
global q = [0.892, 0.919, 0.953, 0.914, 0.709, 0.88, 0.737, 0.593, 0.903, 0.767, 0.635, 0.191, 0.969, 0.86, 0.689, 0.327, 0.94, 0.8, 0.968, 0.987, 0.888, 0.752, 0.622, 0.957, 0.559, 0.831, 0.802, 0.733, 0.503, 0.434, 0.793, 0.897, 0.745, 0.736, 0.401, 0.274, 0.989, 0.998, 0.46, 0.855, 0.804, 0.991, 0.353, 0.688, 0.664, 0.549, 0.803, 0.933, 0.994, 0.818, 0.893, 0.608, 0.933, 0.728, 0.984, 0.79, 0.994, 0.806, 0.845, 0.948, 0.988, 0.815, 0.998, 0.78, 0.901, 0.572, 0.681, 0.205, 0.969, 0.601, 0.983, 0.839, 0.664, 0.993, 0.844, 0.833, 0.772, 0.923, 0.887, 0.733, 0.214, 0.738, 0.856, 0.96, 0.77, 0.797, 0.849, 0.994, 0.828, 0.558, 0.457, 0.876, 0.778, 0.94, 0.998, 0.956, 0.661, 0.792, 0.896, 0.689, 0.97, 0.705, 0.658, 0.598, 0.443, 0.584, 0.895, 0.925, 0.632, 0.951, 0.857, 0.933, 0.049, 0.9, 0.726, 0.98, 0.864, 0.826, 0.763, 0.165, 0.991, 0.733, 0.743, 0.87, 0.692, 0.445, 0.921, 0.982, 0.992, 0.41, 0.921, 0.999, 0.74, 0.741, 0.727, 0.986, 0.998, 0.964, 0.273, 0.785, 0.983, 0.597, 0.793, 0.764, 0.7, 0.942, 0.283, 0.353, 0.789, 0.942, 0.181, 0.763, 0.759, 0.94, 0.9, 0.892, 0.694, 0.987, 0.936, 0.914, 0.494, 0.918, 0.863, 0.498, 0.691, 0.977, 0.72, 0.378, 0.814, 0.651, 0.759, 0.86, 0.905, 0.952, 0.584, 0.338, 0.861, 0.728, 0.911, 0.82, 0.787, 0.98, 0.322, 0.877, 0.915, 0.921, 0.985, 0.852, 0.926, 0.937, 0.72, 0.274, 0.913, 0.898, 0.729, 0.766, 0.928, 0.646, 0.599, 0.993, 0.666, 0.946, 0.948, 0.769, 0.604, 0.848, 0.923, 0.752, 0.853, 0.579, 0.763, 0.657, 0.652, 0.687, 0.999, 0.975, 0.781, 0.782, 0.958, 0.863, 0.52, 0.302, 0.987, 0.457, 0.864]
global origin = 1
global destination = 50