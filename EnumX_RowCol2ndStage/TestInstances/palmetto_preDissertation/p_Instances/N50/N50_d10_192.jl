global arcs = [1 10; 1 34; 1 47; 1 48; 1 49; 1 50; 2 7; 2 13; 2 17; 2 22; 2 35; 2 46; 2 47; 2 50; 3 7; 3 11; 3 18; 3 24; 3 42; 3 44; 3 49; 4 17; 4 18; 4 23; 4 26; 4 34; 4 38; 5 11; 5 14; 5 17; 5 28; 5 31; 5 34; 5 46; 6 5; 6 17; 6 32; 6 33; 6 35; 6 40; 7 5; 7 11; 7 13; 7 17; 7 35; 7 43; 8 9; 8 17; 8 22; 8 23; 8 24; 8 37; 9 7; 9 19; 9 28; 9 42; 10 20; 10 21; 10 23; 10 28; 10 46; 10 48; 10 50; 11 13; 11 17; 11 20; 11 26; 11 46; 12 5; 12 9; 12 11; 12 37; 13 9; 13 17; 13 30; 13 33; 13 39; 13 43; 14 4; 14 20; 14 23; 14 32; 14 37; 14 46; 15 13; 16 2; 16 15; 16 24; 16 26; 16 29; 16 30; 16 34; 16 50; 17 6; 17 33; 17 36; 17 45; 17 48; 18 29; 18 30; 18 34; 18 39; 19 2; 19 12; 19 23; 19 29; 19 37; 20 6; 20 13; 20 17; 20 33; 20 34; 21 5; 21 9; 21 16; 21 20; 21 29; 21 38; 21 47; 22 8; 22 12; 22 16; 22 19; 22 25; 22 50; 23 12; 23 19; 23 20; 23 29; 23 31; 23 35; 23 41; 23 49; 24 30; 24 38; 24 49; 25 2; 25 8; 25 10; 25 27; 25 36; 25 48; 25 49; 26 3; 26 8; 26 9; 26 10; 26 23; 26 25; 26 27; 26 37; 26 38; 26 42; 26 48; 26 49; 27 17; 27 20; 27 30; 27 43; 28 12; 28 19; 28 35; 28 37; 28 38; 29 16; 29 25; 29 32; 29 40; 29 44; 30 9; 30 11; 30 24; 30 33; 30 36; 31 16; 31 18; 31 44; 32 4; 32 12; 32 33; 32 40; 32 43; 33 2; 33 11; 33 48; 34 7; 34 22; 35 3; 35 31; 35 34; 35 41; 35 45; 35 47; 36 6; 36 15; 36 43; 36 45; 37 8; 37 25; 37 27; 37 32; 37 35; 38 14; 38 25; 38 34; 38 39; 38 41; 38 44; 39 2; 39 3; 39 7; 39 13; 39 15; 40 43; 41 9; 41 13; 41 34; 41 45; 42 17; 42 47; 43 3; 43 20; 43 21; 43 28; 43 38; 43 39; 43 45; 44 4; 44 25; 44 36; 44 39; 44 41; 44 45; 45 13; 45 20; 45 49; 46 3; 46 26; 46 36; 46 43; 47 8; 47 29; 48 5; 48 14; 48 21; 48 31; 48 37; 48 46; 49 6; 49 21; 49 38; 49 41; 49 43]
global d_x = [6.0, 8.0, 4.0, 5.0, 3.0, 10.0, 8.0, 5.0, 5.0, 6.0, 7.0, 3.0, 5.0, 10.0, 1.0, 10.0, 4.0, 1.0, 5.0, 10.0, 4.0, 1.0, 5.0, 9.0, 4.0, 5.0, 9.0, 4.0, 5.0, 6.0, 2.0, 3.0, 5.0, 3.0, 5.0, 9.0, 10.0, 2.0, 7.0, 6.0, 7.0, 8.0, 8.0, 3.0, 5.0, 4.0, 9.0, 2.0, 10.0, 6.0, 4.0, 2.0, 6.0, 3.0, 5.0, 10.0, 1.0, 9.0, 3.0, 10.0, 8.0, 5.0, 6.0, 9.0, 10.0, 5.0, 8.0, 7.0, 8.0, 4.0, 1.0, 5.0, 5.0, 5.0, 6.0, 4.0, 2.0, 6.0, 1.0, 6.0, 6.0, 6.0, 8.0, 6.0, 3.0, 5.0, 10.0, 8.0, 10.0, 6.0, 10.0, 6.0, 6.0, 1.0, 5.0, 6.0, 2.0, 2.0, 9.0, 1.0, 5.0, 5.0, 5.0, 9.0, 6.0, 4.0, 3.0, 10.0, 10.0, 4.0, 4.0, 3.0, 3.0, 9.0, 10.0, 8.0, 7.0, 2.0, 4.0, 4.0, 1.0, 7.0, 3.0, 8.0, 4.0, 7.0, 9.0, 10.0, 10.0, 10.0, 3.0, 7.0, 10.0, 7.0, 2.0, 3.0, 5.0, 1.0, 1.0, 10.0, 10.0, 8.0, 8.0, 9.0, 7.0, 8.0, 4.0, 4.0, 2.0, 4.0, 1.0, 7.0, 4.0, 3.0, 5.0, 5.0, 7.0, 5.0, 4.0, 1.0, 10.0, 3.0, 6.0, 7.0, 3.0, 6.0, 9.0, 6.0, 2.0, 8.0, 1.0, 4.0, 4.0, 9.0, 10.0, 1.0, 7.0, 7.0, 5.0, 2.0, 9.0, 1.0, 3.0, 5.0, 2.0, 5.0, 10.0, 2.0, 9.0, 5.0, 10.0, 7.0, 4.0, 9.0, 7.0, 3.0, 3.0, 1.0, 1.0, 7.0, 1.0, 6.0, 3.0, 7.0, 3.0, 8.0, 10.0, 9.0, 9.0, 9.0, 3.0, 8.0, 7.0, 10.0, 9.0, 1.0, 1.0, 7.0, 6.0, 7.0, 2.0, 1.0, 6.0, 9.0, 2.0, 5.0, 5.0, 3.0, 8.0, 6.0, 3.0, 4.0, 7.0, 1.0, 5.0, 3.0, 1.0, 6.0, 5.0, 1.0, 1.0, 7.0, 6.0, 10.0, 10.0, 5.0, 2.0, 7.0, 9.0, 9.0, 4.0, 8.0, 2.0]
global b_x = 5
global d_y = [6.0, 8.0, 6.0, 8.0, 7.0, 4.0, 2.0, 3.0, 6.0, 6.0, 1.0, 1.0, 7.0, 5.0, 5.0, 9.0, 7.0, 3.0, 3.0, 4.0, 4.0, 5.0, 1.0, 8.0, 8.0, 7.0, 5.0, 7.0, 8.0, 3.0, 4.0, 5.0, 9.0, 6.0, 7.0, 7.0, 2.0, 2.0, 7.0, 1.0, 6.0, 3.0, 7.0, 1.0, 3.0, 2.0, 6.0, 1.0, 9.0, 7.0, 5.0, 9.0, 1.0, 9.0, 9.0, 4.0, 4.0, 10.0, 4.0, 4.0, 10.0, 5.0, 1.0, 7.0, 1.0, 4.0, 9.0, 4.0, 3.0, 8.0, 4.0, 10.0, 2.0, 3.0, 1.0, 6.0, 5.0, 10.0, 3.0, 8.0, 8.0, 9.0, 10.0, 8.0, 5.0, 8.0, 6.0, 7.0, 3.0, 2.0, 3.0, 7.0, 10.0, 4.0, 5.0, 6.0, 6.0, 4.0, 4.0, 2.0, 6.0, 6.0, 5.0, 8.0, 1.0, 2.0, 2.0, 9.0, 2.0, 4.0, 7.0, 1.0, 10.0, 2.0, 2.0, 8.0, 4.0, 9.0, 2.0, 5.0, 1.0, 2.0, 2.0, 8.0, 8.0, 9.0, 10.0, 3.0, 4.0, 4.0, 7.0, 3.0, 5.0, 8.0, 3.0, 2.0, 7.0, 9.0, 9.0, 9.0, 9.0, 10.0, 2.0, 8.0, 4.0, 2.0, 6.0, 1.0, 3.0, 3.0, 2.0, 8.0, 7.0, 7.0, 4.0, 1.0, 6.0, 3.0, 9.0, 6.0, 3.0, 6.0, 5.0, 9.0, 1.0, 3.0, 4.0, 5.0, 5.0, 2.0, 5.0, 6.0, 9.0, 3.0, 4.0, 6.0, 2.0, 10.0, 5.0, 2.0, 7.0, 4.0, 4.0, 6.0, 3.0, 5.0, 1.0, 2.0, 1.0, 10.0, 6.0, 8.0, 9.0, 1.0, 4.0, 2.0, 9.0, 5.0, 5.0, 8.0, 8.0, 9.0, 8.0, 4.0, 4.0, 8.0, 10.0, 10.0, 10.0, 5.0, 3.0, 6.0, 7.0, 5.0, 9.0, 4.0, 7.0, 4.0, 9.0, 4.0, 10.0, 5.0, 2.0, 6.0, 2.0, 4.0, 4.0, 2.0, 6.0, 10.0, 6.0, 1.0, 3.0, 8.0, 3.0, 2.0, 3.0, 7.0, 8.0, 9.0, 8.0, 2.0, 8.0, 4.0, 2.0, 2.0, 4.0, 3.0, 10.0, 5.0, 2.0, 10.0, 8.0]
global b_y = 10
global p = [0.464, 0.18, 0.482, 0.041, 0.374, 0.412, 0.572, 0.271, 0.271, 0.402, 0.128, 0.89, 0.644, 0.209, 0.777, 0.065, 0.25, 0.663, 0.203, 0.089, 0.181, 0.478, 0.174, 0.491, 0.721, 0.531, 0.457, 0.099, 0.648, 0.639, 0.135, 0.001, 0.858, 0.431, 0.849, 0.83, 0.698, 0.96, 0.915, 0.534, 0.169, 0.258, 0.603, 0.717, 0.143, 0.685, 0.061, 0.861, 0.639, 0.749, 0.109, 0.884, 0.225, 0.735, 0.739, 0.688, 0.138, 0.778, 0.883, 0.444, 0.461, 0.981, 0.18, 0.898, 0.528, 0.737, 0.934, 0.402, 0.459, 0.199, 0.342, 0.992, 0.862, 0.637, 0.661, 0.663, 0.05, 0.02, 0.531, 0.812, 0.152, 0.969, 0.791, 0.79, 0.18, 0.149, 0.579, 0.09, 0.953, 0.587, 0.637, 0.407, 0.386, 0.478, 0.682, 0.6, 0.398, 0.583, 0.314, 0.711, 0.928, 0.936, 0.49, 0.74, 0.144, 0.805, 0.834, 0.274, 0.439, 0.907, 0.736, 0.738, 0.024, 0.204, 0.664, 0.877, 0.837, 0.401, 0.959, 0.604, 0.871, 0.694, 0.759, 0.347, 0.71, 0.91, 0.039, 0.369, 0.95, 0.264, 0.221, 0.008, 0.894, 0.602, 0.265, 0.917, 0.804, 0.788, 0.071, 0.401, 0.518, 0.031, 0.734, 0.73, 0.132, 0.146, 0.114, 0.236, 0.473, 0.966, 0.196, 0.607, 0.68, 0.736, 0.991, 0.074, 0.772, 0.768, 0.594, 0.867, 0.286, 0.006, 0.681, 0.438, 0.685, 0.748, 0.327, 0.713, 0.341, 0.554, 0.066, 0.689, 0.993, 0.383, 0.199, 0.379, 0.849, 0.637, 0.296, 0.528, 0.956, 0.068, 0.062, 0.848, 0.651, 0.254, 0.305, 0.251, 0.003, 0.142, 0.126, 0.697, 0.747, 0.624, 0.832, 0.635, 0.654, 0.426, 0.907, 0.272, 0.174, 0.68, 0.869, 0.806, 0.781, 0.186, 0.464, 0.818, 0.907, 0.168, 0.852, 0.762, 0.695, 0.481, 0.718, 0.25, 0.794, 0.262, 0.306, 0.706, 0.931, 0.674, 0.661, 0.36, 0.409, 0.421, 0.402, 0.456, 0.545, 0.034, 0.122, 0.947, 0.474, 0.563, 0.785, 0.862, 0.986, 0.768, 0.117, 0.887, 0.235, 0.171, 0.883, 0.638, 0.068, 0.156, 0.057, 0.047, 0.657, 0.76, 0.018, 0.291, 0.306]
global q = [0.868, 0.65, 0.816, 0.411, 0.575, 0.821, 0.767, 0.953, 0.632, 0.9, 0.693, 0.896, 0.982, 0.425, 0.91, 0.464, 0.589, 0.933, 0.647, 0.912, 0.912, 0.935, 0.439, 0.851, 0.762, 0.595, 0.695, 0.86, 0.717, 0.683, 0.255, 0.4, 0.894, 0.527, 0.862, 0.891, 0.832, 0.966, 0.933, 0.624, 0.488, 0.649, 0.616, 0.953, 0.894, 0.894, 0.167, 0.996, 0.783, 0.761, 0.379, 0.978, 0.897, 0.905, 0.769, 0.985, 0.402, 0.842, 0.995, 0.478, 0.909, 0.993, 0.234, 0.967, 0.596, 0.92, 0.962, 0.797, 0.815, 0.715, 0.689, 0.996, 0.962, 0.75, 0.925, 0.902, 0.176, 0.087, 0.628, 0.954, 0.472, 0.975, 0.931, 0.966, 0.672, 0.752, 0.685, 0.493, 0.998, 0.787, 0.923, 0.609, 0.634, 0.95, 0.765, 0.862, 0.679, 0.782, 0.605, 0.99, 0.961, 0.955, 0.889, 0.743, 0.988, 0.969, 0.961, 0.828, 0.512, 0.98, 0.835, 0.997, 0.904, 0.439, 0.724, 0.885, 0.966, 0.92, 0.997, 0.619, 0.959, 0.995, 0.829, 0.681, 0.839, 0.985, 0.445, 0.786, 0.972, 0.941, 0.465, 0.487, 0.901, 0.813, 0.656, 0.92, 0.946, 0.803, 0.263, 0.636, 0.663, 0.476, 0.975, 0.94, 0.939, 0.874, 0.816, 0.312, 0.815, 0.974, 0.495, 0.929, 0.711, 0.77, 0.992, 0.336, 0.885, 0.921, 0.896, 0.951, 0.648, 0.696, 0.971, 0.706, 0.95, 0.965, 0.389, 0.979, 0.455, 0.704, 0.855, 0.853, 0.995, 0.395, 0.64, 0.849, 0.902, 0.824, 0.334, 0.818, 0.969, 0.546, 0.443, 0.971, 0.904, 0.725, 0.6, 0.951, 0.855, 0.191, 0.676, 0.798, 0.984, 0.649, 0.832, 0.834, 0.762, 0.436, 0.994, 0.526, 0.205, 0.724, 0.999, 0.811, 0.821, 0.535, 0.504, 0.955, 0.909, 0.318, 0.976, 0.96, 0.8, 0.775, 0.983, 0.89, 0.896, 0.791, 0.93, 0.866, 0.999, 0.861, 0.681, 0.408, 0.837, 0.732, 0.873, 0.613, 0.892, 0.045, 0.539, 0.988, 0.553, 0.662, 0.878, 0.958, 0.994, 0.844, 0.541, 0.988, 0.598, 0.231, 0.901, 0.816, 0.296, 0.237, 0.441, 0.92, 0.789, 0.874, 0.189, 0.95, 0.464]
global origin = 1
global destination = 50