global arcs = [1 7; 1 16; 1 26; 1 43; 2 6; 2 9; 2 13; 2 31; 2 33; 2 41; 2 49; 3 20; 3 26; 3 42; 4 7; 4 13; 4 23; 4 39; 4 45; 5 9; 5 15; 5 47; 6 40; 7 5; 7 23; 7 31; 7 32; 7 34; 7 42; 7 45; 8 7; 8 23; 8 24; 8 26; 8 40; 9 5; 9 14; 9 23; 9 40; 9 47; 10 5; 10 9; 10 19; 10 21; 10 26; 11 2; 11 14; 11 18; 11 19; 11 20; 11 27; 11 37; 12 3; 12 15; 12 20; 12 42; 12 48; 13 14; 13 15; 13 23; 13 30; 13 36; 13 39; 14 26; 15 9; 15 28; 15 41; 16 19; 16 28; 16 35; 16 47; 17 3; 17 6; 17 7; 17 15; 17 22; 17 25; 17 29; 18 12; 18 27; 18 38; 18 40; 19 20; 19 41; 19 44; 20 32; 20 34; 20 35; 20 40; 21 10; 21 22; 21 28; 21 35; 21 36; 21 46; 22 5; 22 12; 22 32; 22 38; 23 10; 23 13; 23 19; 23 40; 24 5; 24 11; 24 33; 24 48; 24 50; 25 6; 25 14; 25 17; 25 23; 25 33; 25 38; 25 42; 26 7; 27 10; 27 18; 27 26; 27 37; 27 39; 27 45; 28 3; 28 13; 28 18; 28 31; 28 32; 28 44; 29 6; 29 16; 29 20; 29 24; 29 37; 29 39; 30 22; 30 49; 31 2; 31 4; 31 16; 31 23; 32 2; 32 3; 32 9; 32 13; 32 21; 32 37; 32 43; 32 45; 32 49; 33 18; 33 20; 33 23; 33 36; 33 43; 33 48; 34 19; 34 37; 34 50; 35 3; 35 5; 35 9; 35 15; 35 33; 36 3; 36 13; 36 25; 36 28; 36 29; 36 30; 36 34; 36 40; 37 7; 37 19; 37 30; 37 36; 37 45; 37 50; 38 15; 38 21; 38 31; 38 45; 38 50; 39 3; 39 10; 39 15; 39 21; 39 23; 39 45; 40 3; 40 7; 40 9; 40 10; 40 12; 40 29; 40 45; 41 7; 41 9; 41 11; 41 25; 41 28; 41 29; 41 33; 41 34; 41 40; 41 50; 42 20; 42 29; 42 40; 43 20; 43 37; 43 41; 43 44; 43 48; 44 15; 44 23; 44 41; 44 43; 45 6; 45 14; 45 23; 45 32; 46 4; 46 7; 46 8; 46 18; 46 26; 46 27; 46 37; 46 45; 47 13; 47 35; 47 42; 47 48; 48 10; 48 33; 48 35; 48 36; 48 49; 49 9; 49 31; 49 36]
global d_x = [3.0, 6.0, 7.0, 4.0, 5.0, 10.0, 10.0, 8.0, 9.0, 1.0, 1.0, 7.0, 9.0, 8.0, 7.0, 2.0, 4.0, 10.0, 9.0, 6.0, 4.0, 10.0, 6.0, 9.0, 2.0, 10.0, 4.0, 4.0, 1.0, 7.0, 7.0, 4.0, 3.0, 6.0, 8.0, 4.0, 6.0, 5.0, 2.0, 6.0, 7.0, 2.0, 6.0, 1.0, 4.0, 1.0, 10.0, 2.0, 5.0, 8.0, 7.0, 8.0, 3.0, 3.0, 5.0, 9.0, 5.0, 10.0, 5.0, 9.0, 8.0, 3.0, 10.0, 2.0, 1.0, 10.0, 10.0, 1.0, 6.0, 9.0, 5.0, 9.0, 10.0, 10.0, 3.0, 6.0, 6.0, 1.0, 10.0, 2.0, 5.0, 2.0, 7.0, 5.0, 3.0, 2.0, 1.0, 10.0, 10.0, 4.0, 1.0, 7.0, 2.0, 9.0, 3.0, 7.0, 6.0, 3.0, 10.0, 4.0, 8.0, 4.0, 10.0, 10.0, 2.0, 2.0, 7.0, 10.0, 7.0, 7.0, 5.0, 5.0, 9.0, 7.0, 1.0, 4.0, 5.0, 7.0, 2.0, 5.0, 3.0, 4.0, 5.0, 8.0, 2.0, 6.0, 8.0, 3.0, 6.0, 2.0, 2.0, 9.0, 9.0, 6.0, 8.0, 9.0, 5.0, 7.0, 6.0, 7.0, 8.0, 7.0, 7.0, 1.0, 10.0, 5.0, 7.0, 9.0, 7.0, 9.0, 1.0, 2.0, 2.0, 2.0, 1.0, 5.0, 5.0, 6.0, 7.0, 5.0, 1.0, 3.0, 5.0, 7.0, 8.0, 1.0, 3.0, 5.0, 1.0, 6.0, 2.0, 10.0, 5.0, 3.0, 4.0, 3.0, 6.0, 5.0, 2.0, 2.0, 6.0, 10.0, 8.0, 4.0, 6.0, 1.0, 2.0, 4.0, 6.0, 3.0, 6.0, 7.0, 3.0, 8.0, 2.0, 6.0, 10.0, 10.0, 7.0, 6.0, 5.0, 1.0, 3.0, 2.0, 5.0, 10.0, 10.0, 3.0, 6.0, 4.0, 6.0, 7.0, 5.0, 9.0, 3.0, 10.0, 7.0, 3.0, 2.0, 6.0, 4.0, 8.0, 1.0, 3.0, 10.0, 4.0, 7.0, 1.0, 7.0, 10.0, 6.0, 10.0, 2.0, 4.0, 1.0, 1.0, 4.0, 8.0, 2.0, 6.0, 2.0]
global b_x = 5
global d_y = [6.0, 4.0, 3.0, 3.0, 4.0, 9.0, 8.0, 6.0, 2.0, 5.0, 6.0, 6.0, 2.0, 10.0, 2.0, 2.0, 1.0, 1.0, 2.0, 3.0, 8.0, 7.0, 9.0, 9.0, 10.0, 10.0, 4.0, 9.0, 5.0, 7.0, 10.0, 3.0, 3.0, 2.0, 3.0, 6.0, 6.0, 9.0, 8.0, 2.0, 4.0, 3.0, 9.0, 5.0, 9.0, 3.0, 3.0, 7.0, 6.0, 7.0, 10.0, 1.0, 5.0, 8.0, 5.0, 3.0, 9.0, 2.0, 4.0, 10.0, 3.0, 7.0, 5.0, 6.0, 1.0, 9.0, 7.0, 5.0, 9.0, 4.0, 7.0, 3.0, 4.0, 10.0, 10.0, 10.0, 8.0, 9.0, 3.0, 6.0, 5.0, 5.0, 6.0, 2.0, 2.0, 7.0, 9.0, 8.0, 7.0, 3.0, 8.0, 4.0, 8.0, 9.0, 3.0, 7.0, 1.0, 6.0, 4.0, 2.0, 6.0, 7.0, 9.0, 9.0, 10.0, 8.0, 6.0, 9.0, 9.0, 2.0, 4.0, 5.0, 8.0, 1.0, 4.0, 9.0, 7.0, 4.0, 3.0, 4.0, 8.0, 7.0, 10.0, 3.0, 2.0, 7.0, 2.0, 1.0, 8.0, 8.0, 8.0, 3.0, 6.0, 1.0, 10.0, 9.0, 5.0, 3.0, 4.0, 9.0, 5.0, 1.0, 7.0, 10.0, 6.0, 1.0, 3.0, 1.0, 8.0, 5.0, 6.0, 9.0, 7.0, 1.0, 1.0, 5.0, 9.0, 2.0, 6.0, 9.0, 6.0, 1.0, 9.0, 2.0, 5.0, 8.0, 5.0, 10.0, 3.0, 10.0, 10.0, 1.0, 8.0, 4.0, 7.0, 3.0, 4.0, 6.0, 8.0, 9.0, 8.0, 8.0, 10.0, 10.0, 7.0, 10.0, 2.0, 1.0, 5.0, 10.0, 2.0, 3.0, 4.0, 8.0, 2.0, 2.0, 5.0, 7.0, 10.0, 4.0, 4.0, 2.0, 9.0, 1.0, 6.0, 7.0, 4.0, 1.0, 4.0, 2.0, 10.0, 10.0, 6.0, 9.0, 5.0, 6.0, 2.0, 10.0, 2.0, 3.0, 1.0, 9.0, 6.0, 1.0, 3.0, 2.0, 3.0, 2.0, 10.0, 5.0, 1.0, 1.0, 8.0, 7.0, 6.0, 4.0, 9.0, 4.0, 2.0, 8.0, 5.0]
global b_y = 10
global p = [0.479, 0.403, 0.778, 0.51, 0.372, 0.103, 0.295, 0.694, 0.063, 0.357, 0.409, 0.522, 0.349, 0.373, 0.951, 0.528, 0.621, 0.354, 0.719, 0.409, 0.339, 0.416, 0.97, 0.523, 0.92, 0.562, 0.307, 0.312, 0.018, 0.34, 0.476, 0.521, 0.242, 0.272, 0.149, 0.689, 0.917, 0.253, 0.39, 0.467, 0.512, 0.192, 0.699, 0.307, 0.29, 0.087, 0.717, 0.805, 0.909, 0.103, 0.368, 0.63, 0.807, 0.062, 0.165, 0.506, 0.176, 0.794, 0.387, 0.516, 0.23, 0.607, 0.283, 0.861, 0.242, 0.567, 0.177, 0.37, 0.801, 0.177, 0.078, 0.822, 0.733, 0.255, 0.207, 0.71, 0.587, 0.831, 0.359, 0.132, 0.611, 0.212, 0.886, 0.358, 0.761, 0.065, 0.256, 0.854, 0.624, 0.36, 0.685, 0.868, 0.568, 0.12, 0.566, 0.08, 0.585, 0.904, 0.936, 0.67, 0.419, 0.997, 0.752, 0.908, 0.622, 0.807, 0.729, 0.258, 0.601, 0.937, 0.12, 0.159, 0.477, 0.966, 0.707, 0.044, 0.137, 0.948, 0.654, 0.225, 0.929, 0.569, 0.768, 0.176, 0.348, 0.462, 0.619, 0.427, 0.304, 0.303, 0.937, 0.616, 0.684, 0.572, 0.161, 0.122, 0.633, 0.965, 0.376, 0.46, 0.746, 0.843, 0.342, 0.871, 0.075, 0.708, 0.313, 0.526, 0.127, 0.903, 0.21, 0.597, 0.379, 0.363, 0.489, 0.492, 0.53, 0.09, 0.835, 0.141, 0.86, 0.318, 0.387, 0.799, 0.887, 0.477, 0.365, 0.644, 0.181, 0.72, 0.583, 0.253, 0.682, 0.881, 0.052, 0.868, 0.949, 0.286, 0.84, 0.527, 0.856, 0.588, 0.803, 0.915, 0.971, 0.095, 0.128, 0.175, 0.534, 0.009, 0.59, 0.417, 0.825, 0.375, 0.83, 0.795, 0.794, 0.942, 0.568, 0.665, 0.845, 0.906, 0.917, 0.505, 0.335, 0.294, 0.685, 0.097, 0.403, 0.539, 0.851, 0.284, 0.773, 0.444, 0.635, 0.467, 0.994, 0.432, 0.18, 0.578, 0.309, 0.468, 0.959, 0.159, 0.651, 0.987, 0.964, 0.182, 0.895, 0.443, 0.096, 0.958, 0.378, 0.42, 0.647, 0.109, 0.142, 0.753, 0.819, 0.613, 0.277]
global q = [0.799, 0.623, 0.951, 0.688, 0.376, 0.691, 0.582, 0.798, 0.165, 0.654, 0.42, 0.692, 0.527, 0.421, 0.994, 0.641, 0.727, 0.992, 0.889, 0.849, 0.401, 0.694, 0.981, 0.989, 0.979, 0.907, 0.635, 0.536, 0.279, 0.611, 0.636, 0.713, 0.939, 0.939, 0.202, 0.929, 0.986, 0.564, 0.59, 0.999, 0.683, 0.821, 0.825, 0.559, 0.325, 0.992, 0.81, 0.885, 0.983, 0.429, 0.791, 0.732, 0.893, 0.13, 0.223, 0.76, 0.473, 0.993, 0.936, 0.531, 0.374, 0.872, 0.479, 0.939, 0.59, 0.829, 0.372, 0.95, 0.827, 0.969, 0.154, 0.989, 0.904, 0.489, 0.323, 0.987, 0.964, 0.94, 0.474, 0.813, 0.761, 0.26, 0.99, 0.866, 0.844, 0.87, 0.675, 0.885, 0.794, 0.862, 0.752, 0.907, 0.944, 0.52, 0.771, 0.348, 0.846, 0.944, 0.948, 0.822, 0.42, 0.997, 0.939, 0.96, 0.687, 0.96, 0.824, 0.435, 0.989, 0.983, 0.507, 0.521, 0.842, 0.988, 0.799, 0.597, 0.801, 0.997, 0.721, 0.243, 0.968, 0.768, 0.909, 0.371, 0.541, 0.743, 0.81, 0.795, 0.373, 0.696, 0.978, 0.746, 0.875, 0.786, 0.682, 0.233, 0.975, 0.975, 0.685, 0.897, 0.76, 0.901, 0.89, 0.887, 0.163, 0.933, 0.582, 0.634, 0.36, 0.942, 0.463, 0.785, 0.799, 0.474, 0.925, 0.648, 0.692, 0.425, 0.867, 0.264, 0.932, 0.72, 0.626, 0.828, 0.952, 0.775, 0.495, 0.944, 0.761, 0.83, 0.831, 0.599, 0.726, 0.993, 0.264, 0.899, 0.972, 0.511, 0.979, 0.994, 0.889, 0.607, 0.969, 0.978, 0.997, 0.822, 0.514, 0.842, 0.764, 0.88, 0.636, 0.832, 0.88, 0.512, 0.862, 0.873, 0.998, 0.947, 0.865, 0.8, 0.899, 0.919, 0.917, 0.79, 0.824, 0.427, 0.92, 0.154, 0.645, 0.915, 0.872, 0.465, 0.95, 0.8, 0.9, 0.63, 0.994, 0.763, 0.863, 0.684, 0.886, 0.613, 0.968, 0.236, 0.942, 0.994, 0.968, 0.916, 0.971, 0.505, 0.513, 0.976, 0.887, 0.606, 0.758, 0.328, 0.459, 0.756, 0.999, 0.882, 0.71]
global origin = 1
global destination = 50