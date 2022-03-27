global arcs = [1 27; 1 29; 1 41; 2 9; 2 27; 2 30; 2 40; 3 10; 3 41; 3 46; 3 49; 4 5; 4 14; 4 35; 4 40; 5 6; 5 23; 5 27; 5 30; 5 33; 5 38; 5 43; 5 45; 6 4; 6 7; 6 38; 7 2; 7 8; 7 19; 7 23; 7 25; 7 27; 7 31; 7 41; 8 7; 8 19; 8 22; 8 29; 8 31; 8 33; 8 37; 8 48; 9 11; 9 15; 10 5; 10 18; 11 9; 11 18; 11 23; 11 30; 12 9; 12 16; 12 36; 12 39; 12 41; 13 7; 13 16; 13 25; 13 31; 13 32; 13 41; 14 18; 14 29; 14 42; 15 4; 15 10; 15 17; 15 18; 15 27; 15 34; 15 43; 16 11; 16 29; 16 35; 16 46; 16 48; 17 6; 17 10; 17 13; 17 16; 17 19; 17 21; 17 32; 17 40; 18 3; 18 11; 18 17; 18 34; 18 37; 18 39; 18 40; 18 42; 18 43; 18 44; 18 50; 19 2; 19 3; 19 11; 19 24; 19 26; 19 35; 20 13; 20 25; 20 46; 21 2; 21 5; 21 23; 21 35; 21 46; 22 2; 22 10; 22 14; 22 18; 22 19; 22 23; 22 27; 22 36; 22 37; 23 16; 23 17; 23 20; 23 30; 23 43; 23 46; 23 48; 24 21; 24 27; 24 38; 25 20; 25 26; 25 34; 26 2; 26 5; 26 18; 26 22; 26 43; 26 44; 26 47; 27 12; 27 21; 27 23; 27 25; 27 30; 27 41; 28 4; 28 5; 28 7; 28 8; 28 31; 28 38; 28 45; 28 48; 29 3; 29 7; 29 8; 29 21; 29 32; 29 36; 29 43; 29 48; 30 21; 30 22; 30 23; 30 48; 31 6; 31 7; 31 10; 31 14; 31 22; 31 23; 31 36; 31 48; 31 49; 32 5; 32 21; 32 28; 32 33; 32 36; 32 48; 33 2; 33 5; 33 35; 34 4; 34 23; 34 25; 34 33; 35 3; 35 10; 35 16; 35 48; 36 4; 36 6; 36 10; 36 16; 37 21; 37 22; 37 30; 37 36; 37 38; 37 39; 37 45; 38 7; 38 11; 38 41; 38 47; 39 3; 39 13; 39 28; 39 49; 40 11; 40 27; 40 36; 40 43; 41 7; 41 15; 41 21; 41 26; 41 28; 41 29; 41 31; 41 37; 42 8; 42 12; 42 20; 42 21; 42 43; 43 3; 43 9; 43 11; 43 29; 44 5; 44 15; 44 21; 45 5; 45 18; 45 23; 45 28; 45 38; 46 12; 46 24; 46 25; 46 28; 46 49; 47 4; 47 6; 47 8; 47 22; 47 33; 47 38; 47 50; 48 7; 48 16; 48 34; 49 5; 49 7; 49 16; 49 24; 49 25; 49 36; 49 38]
global d_x = [10.0, 2.0, 1.0, 1.0, 5.0, 10.0, 8.0, 6.0, 3.0, 10.0, 8.0, 2.0, 4.0, 5.0, 10.0, 10.0, 5.0, 10.0, 1.0, 1.0, 2.0, 4.0, 8.0, 10.0, 9.0, 2.0, 6.0, 4.0, 4.0, 10.0, 7.0, 7.0, 4.0, 2.0, 3.0, 1.0, 4.0, 4.0, 6.0, 4.0, 9.0, 4.0, 4.0, 9.0, 6.0, 6.0, 9.0, 6.0, 9.0, 8.0, 6.0, 4.0, 2.0, 10.0, 6.0, 1.0, 2.0, 7.0, 1.0, 2.0, 9.0, 2.0, 8.0, 7.0, 5.0, 9.0, 1.0, 1.0, 1.0, 5.0, 9.0, 6.0, 9.0, 7.0, 3.0, 9.0, 1.0, 5.0, 9.0, 9.0, 1.0, 1.0, 8.0, 6.0, 4.0, 6.0, 8.0, 3.0, 10.0, 5.0, 9.0, 6.0, 10.0, 9.0, 4.0, 5.0, 4.0, 6.0, 3.0, 1.0, 4.0, 5.0, 5.0, 7.0, 3.0, 5.0, 1.0, 9.0, 3.0, 7.0, 8.0, 6.0, 1.0, 6.0, 9.0, 9.0, 10.0, 5.0, 10.0, 7.0, 1.0, 8.0, 4.0, 2.0, 6.0, 8.0, 3.0, 10.0, 7.0, 1.0, 2.0, 9.0, 9.0, 7.0, 10.0, 9.0, 3.0, 10.0, 8.0, 6.0, 3.0, 2.0, 4.0, 6.0, 7.0, 7.0, 5.0, 10.0, 2.0, 3.0, 3.0, 8.0, 1.0, 3.0, 9.0, 5.0, 1.0, 8.0, 8.0, 8.0, 8.0, 1.0, 4.0, 4.0, 3.0, 10.0, 4.0, 5.0, 10.0, 2.0, 5.0, 5.0, 1.0, 7.0, 3.0, 10.0, 3.0, 10.0, 8.0, 4.0, 3.0, 1.0, 6.0, 1.0, 6.0, 1.0, 1.0, 2.0, 2.0, 5.0, 7.0, 4.0, 4.0, 1.0, 8.0, 10.0, 1.0, 10.0, 2.0, 2.0, 4.0, 4.0, 6.0, 10.0, 7.0, 7.0, 5.0, 3.0, 4.0, 4.0, 6.0, 7.0, 10.0, 3.0, 3.0, 2.0, 9.0, 1.0, 5.0, 2.0, 6.0, 7.0, 7.0, 10.0, 4.0, 3.0, 4.0, 10.0, 5.0, 5.0, 3.0, 10.0, 5.0, 7.0, 2.0, 2.0, 5.0, 1.0, 6.0, 1.0, 7.0, 8.0, 8.0, 1.0, 7.0, 7.0, 1.0, 5.0, 10.0, 4.0, 4.0, 10.0, 4.0, 5.0, 7.0, 5.0, 3.0, 2.0, 10.0, 1.0]
global b_x = 5
global d_y = [4.0, 3.0, 1.0, 2.0, 4.0, 10.0, 7.0, 8.0, 5.0, 9.0, 8.0, 4.0, 3.0, 10.0, 9.0, 4.0, 7.0, 9.0, 2.0, 1.0, 10.0, 10.0, 6.0, 7.0, 3.0, 6.0, 2.0, 3.0, 9.0, 4.0, 1.0, 7.0, 8.0, 3.0, 1.0, 1.0, 7.0, 9.0, 4.0, 10.0, 5.0, 5.0, 9.0, 4.0, 6.0, 4.0, 2.0, 10.0, 4.0, 4.0, 2.0, 8.0, 8.0, 8.0, 1.0, 2.0, 5.0, 10.0, 7.0, 6.0, 2.0, 5.0, 5.0, 9.0, 10.0, 7.0, 2.0, 9.0, 7.0, 3.0, 3.0, 3.0, 7.0, 3.0, 4.0, 4.0, 5.0, 8.0, 6.0, 1.0, 5.0, 5.0, 3.0, 1.0, 1.0, 6.0, 9.0, 5.0, 7.0, 5.0, 4.0, 6.0, 4.0, 3.0, 2.0, 8.0, 10.0, 4.0, 1.0, 4.0, 1.0, 2.0, 1.0, 3.0, 9.0, 2.0, 10.0, 10.0, 6.0, 8.0, 1.0, 2.0, 6.0, 2.0, 8.0, 5.0, 9.0, 9.0, 5.0, 3.0, 1.0, 10.0, 2.0, 5.0, 6.0, 10.0, 4.0, 2.0, 1.0, 3.0, 10.0, 2.0, 8.0, 1.0, 1.0, 3.0, 6.0, 7.0, 7.0, 9.0, 6.0, 7.0, 10.0, 5.0, 3.0, 8.0, 8.0, 7.0, 10.0, 10.0, 4.0, 1.0, 2.0, 2.0, 1.0, 1.0, 6.0, 9.0, 3.0, 7.0, 10.0, 10.0, 2.0, 2.0, 4.0, 2.0, 9.0, 9.0, 3.0, 3.0, 8.0, 9.0, 4.0, 3.0, 6.0, 10.0, 8.0, 6.0, 5.0, 4.0, 1.0, 10.0, 5.0, 2.0, 1.0, 5.0, 8.0, 2.0, 8.0, 9.0, 10.0, 2.0, 10.0, 7.0, 10.0, 10.0, 3.0, 8.0, 10.0, 4.0, 1.0, 4.0, 10.0, 5.0, 10.0, 1.0, 1.0, 2.0, 6.0, 3.0, 1.0, 7.0, 4.0, 6.0, 1.0, 7.0, 6.0, 6.0, 6.0, 4.0, 7.0, 8.0, 9.0, 10.0, 2.0, 7.0, 3.0, 8.0, 7.0, 9.0, 1.0, 3.0, 4.0, 9.0, 5.0, 1.0, 9.0, 1.0, 8.0, 2.0, 4.0, 10.0, 4.0, 5.0, 10.0, 5.0, 4.0, 5.0, 1.0, 1.0, 10.0, 4.0, 10.0, 3.0, 6.0, 7.0, 3.0, 1.0, 1.0, 1.0]
global b_y = 10
global p = [0.989, 0.121, 0.295, 0.343, 0.859, 0.859, 0.083, 0.405, 0.51, 0.626, 0.062, 0.076, 0.554, 0.922, 0.118, 0.255, 0.459, 0.728, 0.325, 0.11, 0.384, 0.008, 0.791, 0.085, 0.172, 0.044, 0.513, 0.981, 0.015, 0.229, 0.65, 0.166, 0.208, 0.904, 0.839, 0.081, 0.227, 0.252, 0.538, 0.356, 0.929, 0.265, 0.684, 0.746, 0.256, 0.778, 0.786, 0.747, 0.146, 0.41, 0.449, 0.359, 0.7, 0.224, 0.529, 0.418, 0.67, 0.592, 0.084, 0.21, 0.113, 0.234, 0.27, 0.559, 0.677, 0.406, 0.917, 0.371, 0.832, 0.632, 0.717, 0.204, 0.207, 0.01, 0.689, 0.454, 0.137, 0.703, 0.341, 0.77, 0.545, 0.219, 0.735, 0.457, 0.549, 0.266, 0.175, 0.801, 0.639, 0.5, 0.447, 0.829, 0.359, 0.289, 0.29, 0.139, 0.259, 0.563, 0.513, 0.813, 0.279, 0.444, 0.053, 0.094, 0.423, 0.95, 0.792, 0.651, 0.119, 0.451, 0.847, 0.769, 0.038, 0.121, 0.592, 0.709, 0.941, 0.338, 0.621, 0.795, 0.227, 0.807, 0.637, 0.528, 0.495, 0.922, 0.823, 0.491, 0.878, 0.631, 0.576, 0.714, 0.581, 0.328, 0.395, 0.613, 0.003, 0.238, 0.674, 0.803, 0.905, 0.795, 0.622, 0.273, 0.788, 0.274, 0.644, 0.537, 0.756, 0.344, 0.382, 0.512, 0.461, 0.317, 0.129, 0.8, 0.141, 0.85, 0.29, 0.243, 0.912, 0.119, 0.681, 0.533, 0.835, 0.716, 0.605, 0.811, 0.148, 0.044, 0.05, 0.005, 0.153, 0.662, 0.087, 0.572, 0.2, 0.768, 0.519, 0.034, 0.597, 0.619, 0.476, 0.168, 0.036, 0.697, 0.245, 0.391, 0.655, 0.418, 0.043, 0.065, 0.313, 0.371, 0.647, 0.048, 0.064, 0.471, 0.055, 0.252, 0.358, 0.723, 0.81, 0.679, 0.722, 0.13, 0.064, 0.35, 0.312, 0.464, 0.078, 0.586, 0.42, 0.746, 0.126, 0.305, 0.131, 0.648, 0.596, 0.783, 0.327, 0.881, 0.388, 0.509, 0.509, 0.087, 0.742, 0.593, 0.709, 0.55, 0.517, 0.224, 0.821, 0.703, 0.336, 0.238, 0.755, 0.063, 0.431, 0.81, 0.094, 0.833, 0.533, 0.676, 0.682, 0.777, 0.765, 0.61, 0.865, 0.39, 0.056, 0.671, 0.264, 0.603, 0.757, 0.135, 0.334, 0.383, 0.977, 0.151]
global q = [0.998, 0.688, 0.549, 0.741, 0.987, 0.899, 0.509, 0.78, 0.563, 0.923, 0.478, 0.336, 0.623, 0.966, 0.975, 0.505, 0.843, 0.809, 0.609, 0.662, 0.958, 0.898, 0.871, 0.157, 0.788, 0.491, 0.86, 0.989, 0.888, 0.866, 0.926, 0.259, 0.645, 0.933, 0.993, 0.458, 0.734, 0.954, 0.552, 0.576, 0.996, 0.913, 0.777, 0.999, 0.662, 0.913, 0.852, 0.832, 0.694, 0.913, 0.92, 0.77, 0.755, 0.796, 0.987, 0.705, 0.9, 0.896, 0.664, 0.521, 0.547, 0.362, 0.794, 0.995, 0.928, 0.518, 0.95, 0.718, 0.871, 0.676, 0.899, 0.357, 0.257, 0.871, 0.991, 0.636, 0.369, 0.886, 0.933, 0.937, 0.903, 0.887, 0.892, 0.815, 0.761, 0.898, 0.921, 0.939, 0.665, 0.822, 0.922, 0.843, 0.759, 0.974, 0.476, 0.155, 0.27, 0.679, 0.993, 0.868, 0.656, 0.894, 0.789, 0.759, 0.6, 0.972, 0.795, 0.744, 0.337, 0.536, 0.874, 0.997, 0.623, 0.89, 0.698, 0.77, 0.95, 0.425, 0.706, 0.892, 0.895, 0.971, 0.725, 0.657, 0.993, 0.924, 0.902, 0.992, 0.894, 0.907, 0.729, 0.899, 0.768, 0.409, 0.89, 0.823, 0.285, 0.313, 0.787, 0.988, 0.974, 0.946, 0.639, 0.815, 0.876, 0.274, 0.795, 0.905, 0.84, 0.991, 0.627, 0.544, 0.856, 0.492, 0.637, 0.932, 0.773, 0.946, 0.551, 0.549, 0.982, 0.224, 0.688, 0.777, 0.862, 0.941, 0.805, 0.881, 0.697, 0.82, 0.437, 0.441, 0.94, 0.684, 0.147, 0.598, 0.223, 0.999, 0.94, 0.744, 0.863, 0.753, 0.826, 0.837, 0.568, 0.979, 0.857, 0.66, 0.96, 0.902, 0.578, 0.304, 0.431, 0.537, 0.749, 0.181, 0.234, 0.755, 0.672, 0.397, 0.658, 0.827, 0.881, 0.851, 0.836, 0.382, 0.123, 0.812, 0.509, 0.625, 0.223, 0.968, 0.466, 0.946, 0.884, 0.856, 0.923, 0.853, 0.627, 0.854, 0.791, 0.998, 0.83, 0.98, 0.902, 0.125, 0.798, 0.656, 0.872, 0.797, 0.801, 0.865, 0.873, 0.727, 0.443, 0.843, 0.837, 0.303, 0.481, 0.885, 0.899, 0.898, 0.648, 0.812, 0.902, 0.871, 0.986, 0.724, 0.883, 0.634, 0.32, 0.761, 0.579, 0.856, 0.895, 0.526, 0.842, 0.921, 0.991, 0.528]
global origin = 1
global destination = 50