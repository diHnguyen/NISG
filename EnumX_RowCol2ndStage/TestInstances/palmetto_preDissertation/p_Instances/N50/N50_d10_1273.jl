global arcs = [1 15; 1 24; 1 34; 2 28; 2 33; 3 7; 3 47; 4 5; 4 34; 4 35; 5 46; 6 8; 6 11; 6 17; 6 21; 6 35; 6 47; 6 48; 6 50; 7 22; 7 27; 7 34; 7 47; 8 13; 8 32; 8 42; 9 13; 9 19; 9 38; 9 43; 10 18; 10 26; 10 36; 10 40; 10 49; 11 3; 11 16; 11 20; 11 31; 11 39; 11 43; 11 46; 12 29; 13 29; 13 37; 13 41; 14 4; 14 6; 14 27; 14 48; 15 13; 15 19; 15 21; 15 31; 15 37; 16 21; 16 28; 16 32; 16 41; 17 11; 17 19; 17 22; 17 41; 18 2; 18 8; 18 10; 18 27; 18 33; 18 49; 19 22; 19 31; 19 37; 19 41; 19 45; 19 47; 20 2; 20 18; 20 19; 20 39; 20 42; 21 13; 21 36; 22 14; 22 27; 22 46; 23 24; 23 48; 24 6; 24 33; 25 18; 25 32; 25 37; 25 41; 25 44; 26 35; 26 40; 26 49; 27 12; 27 19; 27 28; 27 29; 27 35; 27 38; 28 19; 28 33; 28 34; 28 40; 29 13; 29 31; 29 40; 29 46; 30 9; 30 20; 30 24; 30 35; 30 45; 31 8; 31 12; 31 17; 31 21; 31 29; 31 32; 31 38; 32 6; 32 44; 32 45; 32 49; 33 17; 33 19; 33 23; 33 30; 34 7; 34 15; 34 21; 34 32; 34 33; 34 45; 34 46; 34 48; 35 6; 35 9; 35 25; 35 30; 35 42; 35 47; 35 49; 36 8; 36 11; 36 22; 36 23; 36 30; 36 46; 37 20; 37 26; 37 31; 37 38; 38 4; 38 25; 38 35; 38 37; 39 15; 39 18; 39 34; 39 45; 39 49; 40 4; 40 8; 40 12; 40 34; 40 46; 40 47; 41 20; 41 31; 41 33; 41 39; 41 42; 42 6; 42 10; 42 20; 42 40; 42 50; 43 13; 43 18; 43 26; 43 39; 44 3; 44 20; 44 46; 44 48; 44 50; 45 23; 45 24; 45 27; 45 38; 45 39; 45 41; 46 9; 46 14; 46 41; 47 8; 47 17; 47 18; 47 29; 47 35; 47 38; 48 8; 48 17; 48 18; 48 24; 48 31; 48 43; 48 46; 49 4; 49 14; 49 28]
global d_x = [6.0, 10.0, 6.0, 8.0, 1.0, 8.0, 6.0, 6.0, 1.0, 4.0, 3.0, 2.0, 8.0, 5.0, 5.0, 5.0, 8.0, 10.0, 6.0, 1.0, 6.0, 9.0, 3.0, 6.0, 10.0, 9.0, 1.0, 6.0, 7.0, 5.0, 1.0, 6.0, 3.0, 1.0, 1.0, 10.0, 3.0, 9.0, 2.0, 4.0, 5.0, 3.0, 9.0, 5.0, 7.0, 5.0, 4.0, 2.0, 6.0, 9.0, 6.0, 4.0, 2.0, 5.0, 1.0, 5.0, 2.0, 5.0, 4.0, 7.0, 10.0, 1.0, 3.0, 4.0, 8.0, 3.0, 7.0, 2.0, 3.0, 2.0, 6.0, 5.0, 6.0, 7.0, 7.0, 9.0, 1.0, 6.0, 4.0, 3.0, 6.0, 5.0, 10.0, 6.0, 1.0, 10.0, 8.0, 5.0, 9.0, 6.0, 6.0, 9.0, 8.0, 1.0, 9.0, 5.0, 3.0, 9.0, 8.0, 1.0, 10.0, 6.0, 1.0, 1.0, 2.0, 1.0, 6.0, 8.0, 9.0, 3.0, 2.0, 9.0, 6.0, 1.0, 10.0, 2.0, 7.0, 4.0, 4.0, 5.0, 10.0, 5.0, 4.0, 9.0, 4.0, 10.0, 8.0, 9.0, 4.0, 1.0, 10.0, 8.0, 10.0, 9.0, 1.0, 10.0, 2.0, 10.0, 3.0, 5.0, 9.0, 4.0, 9.0, 8.0, 3.0, 3.0, 6.0, 7.0, 10.0, 5.0, 7.0, 1.0, 9.0, 2.0, 4.0, 10.0, 7.0, 4.0, 2.0, 9.0, 3.0, 9.0, 2.0, 6.0, 9.0, 5.0, 5.0, 10.0, 6.0, 4.0, 1.0, 7.0, 2.0, 6.0, 4.0, 5.0, 2.0, 9.0, 9.0, 2.0, 2.0, 4.0, 5.0, 9.0, 1.0, 6.0, 2.0, 7.0, 8.0, 2.0, 10.0, 7.0, 3.0, 7.0, 4.0, 8.0, 1.0, 4.0, 5.0, 10.0, 5.0, 10.0, 1.0, 7.0, 6.0, 8.0, 8.0, 5.0, 6.0, 7.0, 9.0, 3.0, 8.0, 5.0, 8.0]
global b_x = 5
global d_y = [2.0, 8.0, 5.0, 10.0, 3.0, 10.0, 6.0, 7.0, 2.0, 1.0, 9.0, 8.0, 9.0, 3.0, 7.0, 9.0, 1.0, 1.0, 8.0, 9.0, 7.0, 6.0, 1.0, 2.0, 9.0, 3.0, 9.0, 3.0, 4.0, 3.0, 2.0, 2.0, 10.0, 2.0, 9.0, 1.0, 2.0, 5.0, 7.0, 1.0, 10.0, 1.0, 2.0, 2.0, 3.0, 8.0, 5.0, 3.0, 6.0, 8.0, 4.0, 8.0, 4.0, 9.0, 5.0, 3.0, 4.0, 1.0, 2.0, 7.0, 1.0, 3.0, 4.0, 5.0, 7.0, 7.0, 9.0, 1.0, 4.0, 8.0, 7.0, 8.0, 4.0, 2.0, 8.0, 10.0, 9.0, 9.0, 7.0, 2.0, 9.0, 2.0, 6.0, 4.0, 7.0, 10.0, 5.0, 1.0, 5.0, 10.0, 9.0, 6.0, 5.0, 3.0, 6.0, 7.0, 10.0, 4.0, 9.0, 2.0, 7.0, 6.0, 7.0, 7.0, 4.0, 2.0, 2.0, 8.0, 2.0, 9.0, 1.0, 10.0, 6.0, 7.0, 6.0, 8.0, 1.0, 10.0, 4.0, 4.0, 1.0, 10.0, 5.0, 3.0, 5.0, 3.0, 8.0, 3.0, 10.0, 6.0, 2.0, 2.0, 6.0, 2.0, 1.0, 7.0, 4.0, 4.0, 2.0, 5.0, 7.0, 2.0, 5.0, 6.0, 7.0, 9.0, 10.0, 8.0, 10.0, 2.0, 10.0, 7.0, 9.0, 4.0, 9.0, 3.0, 3.0, 6.0, 2.0, 2.0, 2.0, 3.0, 10.0, 3.0, 4.0, 5.0, 4.0, 6.0, 3.0, 10.0, 1.0, 7.0, 3.0, 4.0, 9.0, 9.0, 2.0, 8.0, 8.0, 5.0, 9.0, 3.0, 6.0, 9.0, 7.0, 6.0, 1.0, 1.0, 1.0, 7.0, 3.0, 4.0, 3.0, 3.0, 2.0, 8.0, 9.0, 3.0, 1.0, 4.0, 7.0, 4.0, 10.0, 5.0, 5.0, 7.0, 4.0, 1.0, 1.0, 7.0, 10.0, 6.0, 8.0, 4.0, 9.0]
global b_y = 10
global p = [0.902, 0.813, 0.492, 0.431, 0.424, 0.013, 0.945, 0.29, 0.332, 0.268, 0.449, 0.333, 0.562, 0.698, 0.177, 0.937, 0.705, 0.08, 0.72, 0.345, 0.482, 0.835, 0.198, 0.892, 0.666, 0.099, 0.001, 0.705, 0.9, 0.006, 0.526, 0.199, 0.011, 0.502, 0.424, 0.38, 0.013, 0.691, 0.214, 0.072, 0.874, 0.478, 0.434, 0.116, 0.274, 0.125, 0.382, 0.04, 0.695, 0.726, 0.58, 0.192, 0.731, 0.818, 0.986, 0.853, 0.764, 0.477, 0.866, 0.645, 0.24, 0.101, 0.594, 0.389, 0.089, 0.401, 0.362, 0.657, 0.781, 0.523, 0.455, 0.266, 0.141, 0.827, 0.146, 0.285, 0.848, 0.832, 0.152, 0.136, 0.171, 0.126, 0.486, 0.481, 0.46, 0.29, 0.88, 0.476, 0.11, 0.248, 0.77, 0.237, 0.746, 0.125, 0.025, 0.657, 0.157, 0.135, 0.402, 0.074, 0.909, 0.772, 0.154, 0.867, 0.064, 0.112, 0.068, 0.953, 0.335, 0.868, 0.2, 0.913, 0.158, 0.731, 0.075, 0.007, 0.887, 0.443, 0.819, 0.013, 0.783, 0.582, 0.224, 0.489, 0.962, 0.136, 0.609, 0.006, 0.345, 0.778, 0.53, 0.321, 0.016, 0.863, 0.075, 0.281, 0.501, 0.799, 0.974, 0.094, 0.659, 0.585, 0.619, 0.102, 0.797, 0.938, 0.898, 0.66, 0.633, 0.709, 0.986, 0.661, 0.631, 0.73, 0.223, 0.833, 0.705, 0.224, 0.916, 0.936, 0.977, 0.533, 0.761, 0.429, 0.985, 0.332, 0.174, 0.93, 0.05, 0.209, 0.306, 0.797, 0.627, 0.543, 0.378, 0.465, 0.671, 0.603, 0.189, 0.793, 0.726, 0.264, 0.337, 0.436, 0.958, 0.438, 0.609, 0.546, 0.056, 0.608, 0.972, 0.406, 0.643, 0.148, 0.647, 0.1, 0.692, 0.333, 0.256, 0.64, 0.89, 0.598, 0.146, 0.003, 0.254, 0.565, 0.999, 0.058, 0.457, 0.43, 0.786, 0.59, 0.856, 0.285, 0.781]
global q = [0.935, 0.877, 0.552, 0.959, 0.769, 0.664, 0.954, 0.323, 0.476, 0.441, 0.964, 0.391, 0.761, 0.739, 0.639, 0.983, 0.966, 0.087, 0.757, 0.625, 0.727, 0.994, 0.266, 0.959, 0.937, 0.258, 0.874, 0.805, 0.932, 0.314, 0.643, 0.601, 0.251, 0.622, 0.91, 0.562, 0.193, 0.801, 0.712, 0.28, 0.956, 0.939, 0.897, 0.259, 0.693, 0.971, 0.679, 0.391, 0.841, 0.861, 0.709, 0.583, 0.997, 0.982, 0.997, 0.964, 0.875, 0.885, 0.976, 0.882, 0.89, 0.442, 0.656, 0.685, 0.331, 0.812, 0.395, 0.872, 0.814, 0.637, 0.953, 0.798, 0.58, 0.88, 0.245, 0.633, 0.97, 0.935, 0.765, 0.171, 0.669, 0.341, 0.625, 0.723, 0.665, 0.782, 0.966, 0.931, 0.627, 0.274, 0.891, 0.959, 0.781, 0.174, 0.513, 0.662, 0.839, 0.638, 0.945, 0.324, 0.914, 0.951, 0.691, 0.981, 0.21, 0.547, 0.372, 0.969, 0.357, 0.985, 0.305, 0.917, 0.817, 0.924, 0.937, 0.888, 0.954, 0.543, 0.824, 0.642, 0.858, 0.642, 0.98, 0.899, 0.992, 0.616, 0.778, 0.833, 0.879, 0.874, 0.891, 0.495, 0.075, 0.936, 0.444, 0.401, 0.792, 0.97, 0.989, 0.316, 0.794, 0.808, 0.672, 0.893, 0.848, 0.944, 0.964, 0.713, 0.744, 0.798, 0.995, 0.943, 0.883, 0.805, 0.386, 0.864, 0.745, 0.966, 0.922, 0.985, 0.999, 0.935, 0.985, 0.642, 0.986, 0.341, 0.677, 0.933, 0.254, 0.321, 0.641, 0.997, 0.887, 0.855, 0.535, 0.911, 0.836, 0.693, 0.818, 0.943, 0.829, 0.642, 0.935, 0.754, 0.966, 0.659, 0.916, 0.673, 0.883, 0.622, 0.977, 0.935, 0.78, 0.547, 0.863, 0.684, 0.721, 0.958, 0.652, 0.717, 0.993, 0.626, 0.968, 0.701, 0.514, 0.875, 0.999, 0.823, 0.876, 0.838, 0.883, 0.716, 0.858, 0.766, 0.899]
global origin = 1
global destination = 50