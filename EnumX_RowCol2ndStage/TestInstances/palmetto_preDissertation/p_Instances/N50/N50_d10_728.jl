global arcs = [1 7; 1 17; 1 22; 1 41; 2 40; 3 13; 3 17; 3 35; 4 18; 4 43; 5 3; 5 7; 5 14; 5 33; 5 36; 6 3; 6 8; 6 24; 6 29; 6 34; 6 45; 7 5; 7 16; 7 17; 7 22; 7 24; 7 30; 7 42; 7 44; 8 4; 8 47; 9 15; 9 27; 9 36; 9 38; 9 39; 9 43; 9 49; 10 5; 10 14; 11 23; 11 36; 11 38; 11 44; 12 17; 12 34; 12 50; 13 7; 13 19; 13 23; 14 3; 14 13; 14 42; 14 44; 14 49; 15 9; 15 14; 15 21; 15 23; 15 24; 15 34; 15 37; 15 40; 15 48; 15 50; 16 13; 16 17; 16 19; 16 31; 16 39; 16 40; 16 44; 16 48; 17 16; 17 23; 17 31; 17 48; 18 9; 18 38; 19 11; 19 22; 19 41; 20 3; 20 14; 20 29; 20 34; 20 36; 20 45; 20 49; 21 18; 21 40; 21 50; 22 11; 22 27; 22 45; 23 18; 23 21; 23 34; 23 49; 24 10; 24 19; 24 27; 24 28; 24 31; 24 34; 24 37; 25 5; 25 11; 25 29; 25 33; 25 43; 25 46; 26 2; 26 3; 26 27; 26 41; 26 44; 26 48; 27 36; 27 42; 27 45; 27 50; 28 5; 29 10; 29 17; 29 18; 29 20; 29 30; 29 33; 29 39; 29 42; 30 5; 30 6; 30 23; 30 25; 30 36; 30 46; 30 48; 31 11; 31 19; 31 21; 31 35; 31 49; 32 15; 32 26; 32 35; 32 38; 32 40; 32 48; 32 50; 33 13; 33 16; 33 34; 34 2; 34 17; 35 2; 35 3; 35 14; 35 48; 36 27; 36 31; 36 44; 36 46; 37 14; 37 30; 37 32; 37 35; 37 46; 38 5; 38 11; 38 16; 38 17; 38 32; 39 7; 39 15; 39 16; 39 48; 40 4; 40 11; 40 14; 40 21; 40 30; 40 50; 41 26; 41 32; 41 33; 41 35; 41 48; 42 4; 42 6; 42 27; 42 43; 42 48; 43 9; 43 29; 43 36; 44 5; 44 6; 44 15; 44 17; 44 25; 44 47; 44 50; 45 8; 45 17; 45 24; 45 41; 45 43; 45 44; 46 34; 46 37; 46 38; 47 17; 47 19; 47 42; 48 12; 48 13; 48 18; 48 19; 48 39; 48 40; 48 46; 49 4; 49 15]
global d_x = [10.0, 9.0, 5.0, 2.0, 3.0, 1.0, 6.0, 4.0, 7.0, 7.0, 2.0, 5.0, 5.0, 2.0, 7.0, 1.0, 4.0, 10.0, 9.0, 6.0, 9.0, 8.0, 1.0, 3.0, 7.0, 5.0, 10.0, 10.0, 6.0, 2.0, 7.0, 4.0, 4.0, 1.0, 10.0, 10.0, 5.0, 7.0, 8.0, 4.0, 6.0, 2.0, 4.0, 2.0, 1.0, 1.0, 2.0, 5.0, 5.0, 6.0, 9.0, 3.0, 4.0, 6.0, 7.0, 6.0, 2.0, 3.0, 10.0, 3.0, 10.0, 1.0, 10.0, 3.0, 4.0, 3.0, 4.0, 2.0, 8.0, 3.0, 3.0, 4.0, 3.0, 4.0, 4.0, 1.0, 10.0, 10.0, 3.0, 7.0, 3.0, 2.0, 3.0, 7.0, 10.0, 10.0, 3.0, 2.0, 10.0, 10.0, 6.0, 8.0, 2.0, 8.0, 8.0, 6.0, 3.0, 6.0, 9.0, 6.0, 8.0, 2.0, 10.0, 5.0, 6.0, 2.0, 10.0, 5.0, 8.0, 5.0, 1.0, 6.0, 6.0, 8.0, 5.0, 5.0, 3.0, 5.0, 2.0, 3.0, 1.0, 5.0, 3.0, 7.0, 2.0, 3.0, 9.0, 3.0, 3.0, 10.0, 9.0, 7.0, 9.0, 7.0, 7.0, 1.0, 9.0, 9.0, 6.0, 1.0, 10.0, 4.0, 7.0, 5.0, 1.0, 10.0, 2.0, 4.0, 3.0, 10.0, 5.0, 10.0, 6.0, 7.0, 2.0, 6.0, 2.0, 9.0, 9.0, 5.0, 1.0, 5.0, 4.0, 2.0, 1.0, 8.0, 7.0, 10.0, 9.0, 7.0, 5.0, 8.0, 4.0, 1.0, 5.0, 3.0, 5.0, 10.0, 10.0, 3.0, 5.0, 1.0, 10.0, 5.0, 8.0, 1.0, 6.0, 4.0, 9.0, 9.0, 10.0, 9.0, 6.0, 9.0, 4.0, 1.0, 6.0, 5.0, 1.0, 8.0, 1.0, 5.0, 5.0, 5.0, 7.0, 4.0, 5.0, 5.0, 4.0, 4.0, 6.0, 6.0, 5.0, 6.0, 6.0, 8.0, 9.0, 6.0, 6.0, 7.0, 7.0, 4.0, 6.0, 4.0]
global b_x = 5
global d_y = [10.0, 4.0, 5.0, 4.0, 6.0, 1.0, 1.0, 9.0, 9.0, 7.0, 1.0, 1.0, 4.0, 9.0, 2.0, 5.0, 2.0, 8.0, 7.0, 6.0, 4.0, 4.0, 1.0, 3.0, 5.0, 4.0, 2.0, 2.0, 9.0, 2.0, 3.0, 3.0, 5.0, 1.0, 8.0, 1.0, 4.0, 1.0, 2.0, 9.0, 3.0, 10.0, 5.0, 9.0, 4.0, 6.0, 6.0, 4.0, 6.0, 4.0, 8.0, 8.0, 2.0, 5.0, 8.0, 9.0, 8.0, 2.0, 9.0, 6.0, 5.0, 1.0, 9.0, 3.0, 5.0, 3.0, 8.0, 3.0, 1.0, 8.0, 9.0, 8.0, 2.0, 7.0, 6.0, 10.0, 9.0, 2.0, 6.0, 6.0, 1.0, 10.0, 4.0, 4.0, 8.0, 5.0, 1.0, 10.0, 5.0, 5.0, 6.0, 2.0, 5.0, 5.0, 5.0, 10.0, 5.0, 3.0, 8.0, 7.0, 5.0, 8.0, 7.0, 7.0, 2.0, 7.0, 2.0, 1.0, 6.0, 4.0, 8.0, 4.0, 5.0, 1.0, 5.0, 1.0, 6.0, 1.0, 8.0, 10.0, 3.0, 4.0, 10.0, 7.0, 10.0, 3.0, 4.0, 7.0, 1.0, 8.0, 8.0, 6.0, 2.0, 1.0, 10.0, 8.0, 7.0, 1.0, 3.0, 1.0, 6.0, 3.0, 8.0, 3.0, 5.0, 4.0, 10.0, 9.0, 1.0, 2.0, 1.0, 8.0, 2.0, 9.0, 6.0, 6.0, 10.0, 3.0, 8.0, 2.0, 3.0, 7.0, 2.0, 10.0, 5.0, 8.0, 5.0, 1.0, 3.0, 4.0, 2.0, 9.0, 7.0, 8.0, 10.0, 4.0, 1.0, 9.0, 1.0, 9.0, 3.0, 10.0, 2.0, 1.0, 6.0, 10.0, 9.0, 2.0, 2.0, 5.0, 2.0, 9.0, 9.0, 9.0, 7.0, 8.0, 1.0, 6.0, 8.0, 7.0, 5.0, 3.0, 7.0, 8.0, 2.0, 3.0, 8.0, 4.0, 2.0, 9.0, 5.0, 7.0, 7.0, 2.0, 2.0, 6.0, 4.0, 10.0, 4.0, 1.0, 8.0, 2.0, 10.0, 5.0]
global b_y = 10
global p = [0.922, 0.057, 0.714, 0.18, 0.298, 0.734, 0.843, 0.131, 0.575, 0.452, 0.453, 0.478, 0.447, 0.608, 0.511, 0.038, 0.794, 0.483, 0.961, 0.626, 0.472, 0.197, 0.635, 0.396, 0.315, 0.894, 0.465, 0.03, 0.664, 0.344, 0.339, 0.02, 0.082, 0.037, 0.059, 0.453, 0.335, 0.95, 0.923, 0.614, 0.484, 0.11, 0.444, 0.849, 0.354, 0.838, 0.435, 0.332, 0.028, 0.87, 0.241, 0.103, 0.899, 0.188, 0.106, 0.698, 0.392, 0.8, 0.805, 0.198, 0.138, 0.071, 0.132, 0.857, 0.681, 0.963, 0.443, 0.762, 0.572, 0.13, 0.315, 0.756, 0.717, 0.273, 0.344, 0.235, 0.145, 0.566, 0.041, 0.761, 0.468, 0.391, 0.32, 0.302, 0.248, 0.762, 0.454, 0.467, 0.27, 0.169, 0.635, 0.139, 0.857, 0.183, 0.048, 0.192, 0.564, 0.368, 0.979, 0.868, 0.379, 0.492, 0.567, 0.998, 0.974, 0.52, 0.256, 0.327, 0.35, 0.889, 0.47, 0.558, 0.849, 0.628, 0.242, 0.51, 0.112, 0.787, 0.205, 0.513, 0.999, 0.816, 0.815, 0.617, 0.202, 0.681, 0.499, 0.513, 0.984, 0.371, 0.754, 0.802, 0.832, 0.824, 0.79, 0.047, 0.948, 0.005, 0.667, 0.135, 0.499, 0.126, 0.239, 0.627, 0.281, 0.748, 0.692, 0.108, 0.812, 0.163, 0.503, 0.516, 0.394, 0.409, 0.943, 0.746, 0.705, 0.299, 0.145, 0.066, 0.654, 0.515, 0.724, 0.989, 0.782, 0.49, 0.515, 0.944, 0.964, 0.995, 0.959, 0.411, 0.117, 0.689, 0.072, 0.601, 0.116, 0.442, 0.109, 0.399, 0.698, 0.368, 0.875, 0.18, 0.185, 0.673, 0.409, 0.264, 0.451, 0.234, 0.875, 0.98, 0.743, 0.695, 0.047, 0.252, 0.134, 0.656, 0.674, 0.577, 0.441, 0.977, 0.623, 0.816, 0.55, 0.749, 0.622, 0.603, 0.437, 0.379, 0.848, 0.447, 0.895, 0.33, 0.539, 0.681, 0.831, 0.279, 0.021, 0.145, 0.243, 0.316, 0.899, 0.384]
global q = [0.943, 0.809, 0.802, 0.393, 0.94, 0.874, 0.998, 0.949, 0.912, 0.463, 0.759, 0.607, 0.559, 0.904, 0.819, 0.553, 0.887, 0.603, 0.976, 0.885, 0.972, 0.84, 0.647, 0.521, 0.892, 0.993, 0.534, 0.177, 0.932, 0.744, 0.84, 0.476, 0.561, 0.56, 0.304, 0.859, 0.919, 0.961, 0.926, 0.683, 0.648, 0.736, 0.584, 0.907, 0.978, 0.85, 0.635, 0.382, 0.377, 0.955, 0.306, 0.772, 0.95, 0.473, 0.548, 0.935, 0.867, 0.963, 0.836, 0.931, 0.676, 0.415, 0.472, 0.931, 0.953, 0.972, 0.597, 0.808, 0.98, 0.34, 0.718, 0.926, 0.829, 0.414, 0.986, 0.296, 0.828, 0.776, 0.975, 0.875, 0.639, 0.998, 0.425, 0.517, 0.367, 0.83, 0.808, 0.757, 0.804, 0.505, 0.665, 0.298, 0.935, 0.325, 0.341, 0.213, 0.69, 0.983, 0.985, 0.982, 0.846, 0.83, 0.605, 0.999, 0.989, 0.82, 0.305, 0.509, 0.939, 0.895, 0.603, 0.855, 0.994, 0.876, 0.711, 0.521, 0.789, 0.978, 0.841, 0.905, 0.999, 0.82, 0.879, 0.972, 0.823, 0.794, 0.86, 0.833, 0.991, 0.437, 0.794, 0.952, 0.859, 0.884, 0.938, 0.34, 0.975, 0.895, 0.765, 0.148, 0.802, 0.744, 0.935, 0.862, 0.621, 0.875, 0.886, 0.794, 0.863, 0.373, 0.582, 0.765, 0.891, 0.985, 0.948, 0.97, 0.963, 0.328, 0.618, 0.277, 0.986, 0.811, 0.728, 0.992, 0.814, 0.773, 0.717, 0.997, 0.982, 0.995, 0.977, 0.97, 0.98, 0.832, 0.5, 0.687, 0.736, 0.841, 0.828, 0.592, 0.72, 0.721, 0.916, 0.359, 0.421, 0.848, 0.601, 0.303, 0.652, 0.302, 0.976, 0.984, 0.772, 0.745, 0.175, 0.695, 0.927, 0.676, 0.878, 0.777, 0.838, 0.983, 0.776, 0.873, 0.778, 0.931, 0.694, 0.827, 0.889, 0.387, 0.979, 0.857, 0.943, 0.582, 0.709, 0.966, 0.944, 0.353, 0.445, 0.586, 0.983, 0.492, 0.973, 0.777]
global origin = 1
global destination = 50