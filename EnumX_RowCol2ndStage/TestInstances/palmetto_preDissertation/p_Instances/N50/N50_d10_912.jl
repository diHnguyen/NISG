global arcs = [1 23; 1 32; 1 49; 2 14; 2 22; 2 41; 3 4; 3 5; 3 12; 3 16; 3 30; 3 45; 4 15; 4 23; 4 34; 4 35; 4 39; 4 49; 5 3; 5 35; 5 36; 5 37; 5 40; 5 48; 5 49; 6 3; 6 25; 6 34; 6 41; 6 45; 6 47; 6 48; 7 25; 7 33; 8 6; 8 9; 8 16; 8 35; 8 40; 8 43; 9 16; 9 23; 9 26; 9 29; 9 36; 9 37; 9 45; 10 6; 10 24; 10 42; 11 17; 11 25; 11 31; 11 36; 12 7; 12 10; 12 14; 12 17; 12 39; 13 4; 13 22; 13 25; 13 45; 13 47; 14 16; 14 20; 14 24; 14 28; 14 36; 14 37; 15 42; 15 47; 15 50; 16 2; 16 13; 16 27; 16 30; 16 36; 16 38; 16 45; 16 48; 17 5; 17 6; 17 7; 17 27; 17 34; 17 40; 17 47; 18 37; 18 48; 18 50; 19 20; 19 46; 19 50; 20 2; 20 9; 20 24; 20 31; 20 32; 20 36; 20 44; 20 45; 21 38; 21 46; 21 47; 22 6; 22 14; 22 16; 22 38; 23 8; 23 25; 23 47; 23 49; 24 5; 24 10; 24 14; 24 17; 24 25; 24 27; 24 29; 24 44; 24 46; 25 10; 25 13; 25 36; 25 38; 25 45; 25 47; 26 2; 26 8; 26 16; 26 35; 26 37; 26 45; 26 46; 27 35; 27 42; 28 21; 28 38; 28 42; 28 45; 29 9; 29 19; 29 26; 29 36; 29 48; 30 5; 30 27; 30 38; 30 42; 31 26; 31 29; 32 7; 32 13; 32 14; 32 28; 32 36; 32 38; 33 7; 33 24; 33 36; 33 38; 33 50; 34 10; 34 18; 34 21; 34 22; 34 32; 34 35; 34 42; 34 43; 35 21; 35 32; 35 34; 35 42; 35 46; 36 2; 36 3; 36 8; 36 9; 36 13; 36 38; 36 43; 36 46; 37 15; 37 21; 37 42; 37 43; 38 6; 38 20; 38 34; 38 50; 39 26; 39 36; 40 23; 40 37; 40 38; 41 6; 41 12; 41 36; 42 4; 42 9; 42 11; 42 17; 42 24; 42 27; 42 33; 43 24; 43 31; 43 38; 43 46; 44 7; 44 11; 44 14; 44 16; 44 20; 44 26; 44 30; 44 43; 45 2; 45 23; 45 38; 45 43; 45 46; 46 2; 46 3; 46 14; 46 18; 46 28; 46 31; 46 38; 46 44; 46 48; 47 18; 47 29; 47 38; 48 11; 48 42; 49 31]
global d_x = [10.0, 2.0, 2.0, 10.0, 8.0, 2.0, 9.0, 9.0, 4.0, 1.0, 10.0, 4.0, 3.0, 9.0, 7.0, 5.0, 8.0, 3.0, 10.0, 5.0, 10.0, 2.0, 7.0, 1.0, 5.0, 4.0, 2.0, 7.0, 4.0, 8.0, 10.0, 9.0, 8.0, 10.0, 3.0, 3.0, 6.0, 8.0, 4.0, 5.0, 2.0, 1.0, 7.0, 6.0, 5.0, 7.0, 1.0, 3.0, 9.0, 3.0, 8.0, 7.0, 1.0, 5.0, 5.0, 1.0, 9.0, 3.0, 7.0, 5.0, 5.0, 10.0, 1.0, 9.0, 9.0, 4.0, 1.0, 10.0, 4.0, 10.0, 4.0, 10.0, 9.0, 2.0, 10.0, 1.0, 4.0, 9.0, 1.0, 8.0, 6.0, 1.0, 7.0, 4.0, 9.0, 7.0, 6.0, 1.0, 6.0, 8.0, 6.0, 6.0, 8.0, 6.0, 1.0, 10.0, 1.0, 2.0, 9.0, 6.0, 4.0, 4.0, 2.0, 6.0, 9.0, 9.0, 5.0, 7.0, 8.0, 5.0, 10.0, 1.0, 1.0, 1.0, 7.0, 10.0, 10.0, 10.0, 10.0, 4.0, 9.0, 9.0, 2.0, 7.0, 3.0, 9.0, 8.0, 4.0, 3.0, 2.0, 5.0, 4.0, 8.0, 7.0, 7.0, 7.0, 2.0, 10.0, 5.0, 1.0, 3.0, 2.0, 3.0, 1.0, 4.0, 2.0, 9.0, 9.0, 5.0, 6.0, 2.0, 1.0, 6.0, 7.0, 6.0, 7.0, 9.0, 10.0, 6.0, 8.0, 6.0, 9.0, 1.0, 10.0, 8.0, 5.0, 7.0, 6.0, 1.0, 8.0, 7.0, 3.0, 9.0, 1.0, 9.0, 7.0, 7.0, 8.0, 2.0, 1.0, 8.0, 5.0, 6.0, 10.0, 10.0, 8.0, 9.0, 10.0, 1.0, 10.0, 9.0, 1.0, 10.0, 2.0, 4.0, 2.0, 9.0, 5.0, 6.0, 2.0, 10.0, 5.0, 5.0, 3.0, 5.0, 10.0, 8.0, 6.0, 9.0, 1.0, 9.0, 4.0, 8.0, 8.0, 4.0, 4.0, 1.0, 1.0, 9.0, 10.0, 3.0, 9.0, 6.0, 3.0, 9.0, 7.0, 6.0, 8.0, 7.0, 9.0, 8.0, 3.0, 9.0, 7.0, 1.0, 7.0, 7.0, 4.0, 8.0]
global b_x = 5
global d_y = [4.0, 8.0, 1.0, 3.0, 7.0, 6.0, 1.0, 4.0, 9.0, 1.0, 9.0, 2.0, 1.0, 6.0, 6.0, 8.0, 10.0, 3.0, 7.0, 3.0, 3.0, 1.0, 9.0, 10.0, 2.0, 7.0, 8.0, 9.0, 4.0, 6.0, 6.0, 8.0, 8.0, 2.0, 4.0, 10.0, 1.0, 6.0, 2.0, 2.0, 5.0, 6.0, 5.0, 4.0, 4.0, 4.0, 8.0, 10.0, 4.0, 6.0, 9.0, 1.0, 10.0, 1.0, 2.0, 8.0, 6.0, 3.0, 6.0, 8.0, 1.0, 9.0, 5.0, 6.0, 4.0, 8.0, 3.0, 7.0, 8.0, 9.0, 1.0, 10.0, 1.0, 10.0, 9.0, 9.0, 3.0, 1.0, 10.0, 5.0, 4.0, 2.0, 4.0, 4.0, 6.0, 3.0, 5.0, 1.0, 9.0, 9.0, 6.0, 3.0, 1.0, 3.0, 3.0, 6.0, 5.0, 3.0, 2.0, 4.0, 8.0, 4.0, 10.0, 4.0, 6.0, 5.0, 3.0, 2.0, 4.0, 1.0, 5.0, 9.0, 10.0, 1.0, 7.0, 10.0, 6.0, 1.0, 3.0, 4.0, 4.0, 8.0, 5.0, 2.0, 1.0, 3.0, 8.0, 8.0, 7.0, 6.0, 9.0, 3.0, 5.0, 3.0, 7.0, 1.0, 10.0, 2.0, 6.0, 2.0, 2.0, 5.0, 5.0, 1.0, 1.0, 8.0, 1.0, 7.0, 5.0, 4.0, 8.0, 6.0, 4.0, 7.0, 7.0, 1.0, 1.0, 1.0, 3.0, 10.0, 6.0, 3.0, 8.0, 8.0, 9.0, 8.0, 7.0, 5.0, 10.0, 7.0, 9.0, 4.0, 6.0, 8.0, 3.0, 2.0, 8.0, 3.0, 6.0, 7.0, 8.0, 8.0, 1.0, 8.0, 5.0, 7.0, 4.0, 9.0, 2.0, 7.0, 10.0, 8.0, 2.0, 5.0, 5.0, 9.0, 2.0, 1.0, 7.0, 1.0, 7.0, 3.0, 10.0, 3.0, 3.0, 4.0, 1.0, 9.0, 3.0, 5.0, 7.0, 4.0, 1.0, 2.0, 9.0, 6.0, 1.0, 7.0, 8.0, 10.0, 4.0, 3.0, 7.0, 10.0, 6.0, 7.0, 9.0, 4.0, 8.0, 8.0, 2.0, 1.0, 5.0, 8.0, 9.0, 4.0, 6.0, 6.0, 9.0]
global b_y = 10
global p = [0.088, 0.692, 0.101, 0.621, 0.273, 0.887, 0.923, 0.241, 0.84, 0.924, 0.772, 0.648, 0.562, 0.437, 0.055, 0.311, 0.9, 0.178, 0.369, 0.9, 0.91, 0.358, 0.382, 0.921, 0.679, 0.557, 0.749, 0.47, 0.559, 0.211, 0.992, 0.606, 0.029, 0.048, 0.48, 0.163, 0.152, 0.605, 0.724, 0.127, 0.09, 0.355, 0.134, 0.786, 0.566, 0.944, 0.093, 0.961, 0.05, 0.145, 0.885, 0.41, 0.228, 0.783, 0.973, 0.623, 0.129, 0.004, 0.004, 0.452, 0.789, 0.824, 0.4, 0.895, 0.545, 0.685, 0.624, 0.086, 0.719, 0.322, 0.022, 0.8, 0.366, 0.117, 0.479, 0.074, 0.866, 0.981, 0.029, 0.324, 0.397, 0.483, 0.903, 0.106, 0.358, 0.456, 0.91, 0.057, 0.256, 0.404, 0.712, 0.422, 0.615, 0.747, 0.929, 0.725, 0.988, 0.958, 0.785, 0.99, 0.137, 0.926, 0.638, 0.693, 0.717, 0.711, 0.011, 0.637, 0.409, 0.854, 0.173, 0.498, 0.823, 0.7, 0.706, 0.522, 0.822, 0.736, 0.636, 0.222, 0.001, 0.726, 0.841, 0.389, 0.163, 0.513, 0.907, 0.76, 0.297, 0.731, 0.489, 0.197, 0.118, 0.664, 0.517, 0.814, 0.472, 0.362, 0.527, 0.676, 0.829, 0.032, 0.516, 0.728, 0.954, 0.382, 0.251, 0.232, 0.624, 0.224, 0.172, 0.17, 0.651, 0.744, 0.85, 0.977, 0.368, 0.083, 0.077, 0.164, 0.617, 0.583, 0.344, 0.111, 0.198, 0.027, 0.853, 0.536, 0.726, 0.899, 0.83, 0.621, 0.916, 0.695, 0.858, 0.098, 0.934, 0.68, 0.034, 0.106, 0.681, 0.958, 0.251, 0.277, 0.629, 0.958, 0.069, 0.838, 0.42, 0.405, 0.435, 0.178, 0.049, 0.078, 0.893, 0.703, 0.711, 0.168, 0.362, 0.651, 0.131, 0.074, 0.249, 0.714, 0.954, 0.67, 0.939, 0.728, 0.789, 0.083, 0.379, 0.963, 0.231, 0.925, 0.851, 0.527, 0.04, 0.435, 0.601, 0.842, 0.658, 0.267, 0.203, 0.843, 0.836, 0.371, 0.147, 0.78, 0.435, 0.276, 0.557, 0.272, 0.453, 0.408, 0.415, 0.494, 0.851, 0.049, 0.61]
global q = [0.165, 0.97, 0.6, 0.831, 0.633, 0.969, 0.94, 0.541, 0.949, 0.956, 0.835, 0.802, 0.817, 0.687, 0.734, 0.893, 0.982, 0.329, 0.625, 0.947, 0.911, 0.559, 0.907, 0.978, 0.688, 0.778, 0.942, 0.551, 0.755, 0.746, 0.992, 0.987, 0.612, 0.215, 0.731, 0.62, 0.159, 0.816, 0.82, 0.973, 0.806, 0.524, 0.715, 0.922, 0.586, 0.976, 0.572, 0.978, 0.077, 0.499, 0.931, 0.794, 0.279, 0.959, 0.982, 0.864, 0.712, 0.22, 0.086, 0.7, 0.85, 0.962, 0.672, 0.987, 0.563, 0.722, 0.863, 0.217, 0.837, 0.455, 0.956, 0.802, 0.483, 0.145, 0.481, 0.792, 0.924, 0.997, 0.138, 0.799, 0.42, 0.742, 0.926, 0.794, 0.952, 0.525, 0.989, 0.525, 0.7, 0.744, 0.903, 0.476, 0.882, 0.874, 0.931, 0.992, 0.99, 0.985, 0.873, 0.998, 0.213, 0.937, 0.729, 0.804, 0.812, 0.889, 0.416, 0.88, 0.453, 0.962, 0.353, 0.921, 0.92, 0.839, 0.788, 0.581, 0.883, 0.889, 0.896, 0.545, 0.538, 0.897, 0.919, 0.502, 0.985, 0.643, 0.947, 0.866, 0.958, 0.869, 0.533, 0.401, 0.222, 0.679, 0.755, 0.859, 0.513, 0.957, 0.604, 0.873, 0.914, 0.824, 0.823, 0.743, 0.993, 0.683, 0.821, 0.765, 0.863, 0.648, 0.911, 0.653, 0.932, 0.873, 0.946, 0.999, 0.935, 0.576, 0.176, 0.737, 0.828, 0.882, 0.711, 0.81, 0.765, 0.671, 0.997, 0.599, 0.952, 0.902, 0.847, 0.686, 0.927, 0.948, 0.918, 0.957, 0.96, 0.755, 0.795, 0.245, 0.787, 0.964, 0.487, 0.675, 0.815, 0.967, 0.47, 0.851, 0.511, 0.563, 0.762, 0.499, 0.58, 0.91, 0.911, 0.843, 0.982, 0.86, 0.606, 0.964, 0.242, 0.391, 0.527, 0.966, 0.969, 0.725, 0.965, 0.871, 0.856, 0.507, 0.59, 0.968, 0.401, 0.996, 0.872, 0.92, 0.455, 0.68, 0.602, 0.856, 0.998, 0.483, 0.398, 0.848, 0.959, 0.822, 0.154, 0.813, 0.74, 0.427, 0.867, 0.371, 0.512, 0.612, 0.983, 0.564, 0.968, 0.928, 0.716]
global origin = 1
global destination = 50