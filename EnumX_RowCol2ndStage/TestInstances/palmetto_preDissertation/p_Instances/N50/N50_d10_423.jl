global arcs = [1 2; 1 4; 1 10; 1 12; 1 14; 1 21; 1 32; 1 33; 1 44; 2 5; 2 14; 2 22; 2 24; 2 46; 3 14; 3 37; 3 49; 4 8; 4 10; 4 15; 4 17; 4 25; 4 35; 4 49; 4 50; 5 18; 5 21; 5 36; 5 38; 6 9; 6 34; 6 40; 6 43; 6 49; 6 50; 7 9; 7 12; 7 17; 7 19; 7 27; 7 28; 7 42; 7 47; 8 6; 8 13; 8 35; 8 41; 8 48; 8 49; 9 3; 9 12; 9 26; 9 33; 9 42; 9 44; 9 46; 10 2; 10 6; 10 31; 10 49; 11 17; 11 40; 12 6; 12 19; 12 37; 12 49; 13 5; 13 7; 13 17; 13 29; 13 38; 13 48; 14 15; 14 21; 14 22; 14 24; 14 27; 14 30; 14 32; 14 37; 14 41; 14 42; 15 3; 15 16; 15 32; 15 36; 15 37; 15 41; 15 47; 16 4; 16 25; 16 42; 16 45; 17 5; 17 7; 17 23; 17 36; 17 42; 18 36; 19 11; 19 26; 19 47; 20 5; 20 11; 20 21; 20 34; 20 36; 21 5; 21 11; 21 13; 21 15; 21 22; 21 39; 22 3; 22 9; 22 11; 22 13; 22 14; 22 32; 22 37; 23 4; 23 7; 23 42; 24 27; 24 30; 24 34; 24 41; 24 47; 25 3; 25 4; 25 6; 25 19; 25 31; 25 34; 25 40; 26 7; 26 21; 26 50; 27 3; 27 8; 27 22; 27 24; 27 28; 27 45; 28 4; 28 6; 28 7; 28 14; 28 20; 29 6; 29 16; 29 18; 29 19; 30 4; 30 10; 30 18; 30 21; 30 32; 30 46; 31 9; 31 12; 31 13; 31 17; 31 23; 31 30; 31 37; 31 46; 31 48; 32 4; 32 6; 32 20; 32 27; 32 34; 32 48; 33 2; 33 4; 33 8; 33 18; 33 24; 33 47; 34 8; 34 19; 34 21; 34 23; 34 30; 34 35; 34 45; 34 49; 35 12; 35 14; 35 28; 35 40; 35 42; 35 50; 36 31; 37 10; 37 14; 37 33; 37 44; 38 10; 38 17; 38 20; 38 21; 38 39; 38 50; 39 5; 39 15; 39 20; 39 27; 39 45; 39 47; 40 2; 40 11; 40 44; 41 21; 41 25; 41 26; 41 28; 41 31; 41 42; 41 50; 42 15; 42 17; 42 19; 42 22; 42 25; 43 6; 43 19; 43 26; 43 30; 43 36; 44 12; 44 24; 44 29; 44 33; 44 40; 44 50; 45 17; 45 19; 45 20; 45 30; 46 8; 46 27; 46 42; 46 45; 47 8; 47 16; 47 20; 48 4; 48 5; 48 14; 48 28; 48 30; 48 38; 49 25; 49 26; 49 45]
global d_x = [7.0, 4.0, 8.0, 7.0, 8.0, 1.0, 5.0, 8.0, 8.0, 4.0, 2.0, 10.0, 6.0, 10.0, 2.0, 8.0, 8.0, 4.0, 3.0, 9.0, 5.0, 5.0, 6.0, 5.0, 5.0, 3.0, 10.0, 9.0, 8.0, 6.0, 9.0, 5.0, 4.0, 9.0, 7.0, 3.0, 8.0, 6.0, 10.0, 5.0, 1.0, 9.0, 1.0, 6.0, 9.0, 9.0, 2.0, 8.0, 8.0, 2.0, 3.0, 3.0, 10.0, 8.0, 5.0, 6.0, 3.0, 6.0, 1.0, 5.0, 1.0, 6.0, 1.0, 6.0, 1.0, 10.0, 1.0, 9.0, 9.0, 7.0, 1.0, 4.0, 2.0, 6.0, 8.0, 4.0, 8.0, 2.0, 8.0, 1.0, 10.0, 3.0, 5.0, 4.0, 5.0, 3.0, 2.0, 3.0, 5.0, 1.0, 7.0, 7.0, 4.0, 7.0, 7.0, 8.0, 9.0, 4.0, 9.0, 6.0, 10.0, 8.0, 8.0, 4.0, 1.0, 9.0, 6.0, 9.0, 9.0, 2.0, 3.0, 5.0, 4.0, 5.0, 8.0, 1.0, 6.0, 4.0, 9.0, 9.0, 9.0, 4.0, 10.0, 2.0, 10.0, 9.0, 9.0, 2.0, 1.0, 8.0, 10.0, 8.0, 5.0, 2.0, 8.0, 4.0, 3.0, 1.0, 8.0, 3.0, 7.0, 4.0, 1.0, 9.0, 1.0, 8.0, 5.0, 7.0, 7.0, 10.0, 7.0, 4.0, 8.0, 7.0, 3.0, 5.0, 5.0, 2.0, 6.0, 1.0, 2.0, 3.0, 4.0, 10.0, 9.0, 8.0, 8.0, 5.0, 2.0, 10.0, 1.0, 2.0, 10.0, 3.0, 9.0, 3.0, 4.0, 2.0, 5.0, 6.0, 5.0, 2.0, 3.0, 5.0, 5.0, 3.0, 2.0, 9.0, 1.0, 7.0, 6.0, 5.0, 8.0, 9.0, 3.0, 1.0, 2.0, 3.0, 7.0, 3.0, 10.0, 10.0, 5.0, 10.0, 7.0, 3.0, 3.0, 2.0, 7.0, 4.0, 10.0, 8.0, 5.0, 2.0, 3.0, 5.0, 3.0, 1.0, 8.0, 5.0, 4.0, 10.0, 10.0, 6.0, 5.0, 2.0, 5.0, 1.0, 3.0, 7.0, 6.0, 4.0, 5.0, 6.0, 2.0, 8.0, 3.0, 9.0, 6.0, 10.0, 5.0, 9.0, 3.0, 5.0, 6.0, 6.0, 10.0, 8.0, 3.0, 2.0, 6.0, 3.0, 5.0, 1.0, 2.0, 1.0, 1.0]
global b_x = 5
global d_y = [5.0, 6.0, 6.0, 6.0, 5.0, 7.0, 2.0, 8.0, 8.0, 6.0, 1.0, 7.0, 4.0, 10.0, 6.0, 10.0, 3.0, 6.0, 2.0, 2.0, 6.0, 3.0, 8.0, 5.0, 6.0, 2.0, 4.0, 4.0, 7.0, 9.0, 10.0, 1.0, 6.0, 6.0, 8.0, 10.0, 8.0, 10.0, 10.0, 3.0, 3.0, 2.0, 7.0, 9.0, 6.0, 4.0, 6.0, 10.0, 1.0, 6.0, 9.0, 1.0, 5.0, 4.0, 8.0, 7.0, 5.0, 4.0, 3.0, 7.0, 8.0, 4.0, 5.0, 3.0, 1.0, 8.0, 4.0, 4.0, 3.0, 7.0, 4.0, 5.0, 1.0, 5.0, 2.0, 9.0, 3.0, 2.0, 1.0, 10.0, 10.0, 5.0, 7.0, 6.0, 5.0, 10.0, 6.0, 8.0, 7.0, 6.0, 10.0, 6.0, 6.0, 4.0, 2.0, 4.0, 7.0, 1.0, 2.0, 3.0, 6.0, 5.0, 9.0, 1.0, 3.0, 4.0, 4.0, 7.0, 3.0, 10.0, 6.0, 5.0, 8.0, 2.0, 9.0, 7.0, 7.0, 1.0, 7.0, 10.0, 4.0, 4.0, 3.0, 1.0, 2.0, 7.0, 6.0, 8.0, 7.0, 3.0, 2.0, 7.0, 2.0, 3.0, 8.0, 2.0, 7.0, 9.0, 9.0, 4.0, 7.0, 3.0, 3.0, 3.0, 9.0, 10.0, 6.0, 6.0, 3.0, 1.0, 6.0, 6.0, 8.0, 9.0, 5.0, 1.0, 9.0, 6.0, 6.0, 10.0, 9.0, 3.0, 7.0, 9.0, 6.0, 2.0, 4.0, 9.0, 1.0, 4.0, 9.0, 10.0, 7.0, 4.0, 7.0, 4.0, 6.0, 2.0, 7.0, 7.0, 8.0, 4.0, 4.0, 10.0, 9.0, 6.0, 10.0, 9.0, 9.0, 8.0, 1.0, 6.0, 2.0, 7.0, 1.0, 7.0, 4.0, 3.0, 7.0, 8.0, 6.0, 7.0, 4.0, 6.0, 8.0, 9.0, 5.0, 6.0, 2.0, 10.0, 10.0, 4.0, 4.0, 4.0, 9.0, 6.0, 4.0, 6.0, 3.0, 8.0, 7.0, 6.0, 6.0, 2.0, 4.0, 5.0, 10.0, 2.0, 8.0, 8.0, 7.0, 2.0, 8.0, 3.0, 10.0, 10.0, 10.0, 9.0, 8.0, 8.0, 4.0, 6.0, 6.0, 10.0, 8.0, 8.0, 4.0, 6.0, 8.0, 6.0, 4.0, 7.0, 6.0, 5.0, 2.0, 6.0, 7.0]
global b_y = 10
global p = [0.505, 0.632, 0.608, 0.941, 0.416, 0.643, 0.322, 0.84, 0.62, 0.858, 0.531, 0.862, 0.514, 0.706, 0.75, 0.357, 0.138, 0.45, 0.078, 0.595, 0.823, 0.795, 0.396, 0.359, 0.839, 0.547, 0.908, 0.178, 0.688, 0.67, 0.45, 0.123, 0.558, 0.336, 0.907, 0.21, 0.034, 0.979, 0.693, 0.733, 0.161, 0.367, 0.383, 0.265, 0.181, 0.033, 0.286, 0.063, 0.615, 0.865, 0.93, 0.143, 0.562, 0.573, 0.069, 0.971, 0.366, 0.839, 0.77, 0.874, 0.476, 0.361, 0.393, 0.04, 0.137, 0.276, 0.839, 0.839, 0.278, 0.496, 0.977, 0.258, 0.734, 0.567, 0.684, 0.644, 0.569, 0.642, 0.393, 0.841, 0.253, 0.841, 0.702, 0.741, 0.479, 0.592, 0.304, 0.719, 0.981, 0.513, 0.681, 0.388, 0.435, 0.357, 0.396, 0.062, 0.732, 0.736, 0.053, 0.108, 0.685, 0.579, 0.273, 0.977, 0.479, 0.171, 0.92, 0.537, 0.876, 0.661, 0.486, 0.147, 0.295, 0.366, 0.621, 0.905, 0.889, 0.756, 0.798, 0.196, 0.041, 0.508, 0.859, 0.942, 0.086, 0.923, 0.559, 0.546, 0.986, 0.006, 0.203, 0.141, 0.874, 0.768, 0.198, 0.856, 0.777, 0.236, 0.164, 0.02, 0.255, 0.968, 0.949, 0.968, 0.888, 0.389, 0.401, 0.072, 0.955, 0.614, 0.712, 0.753, 0.246, 0.742, 0.219, 0.814, 0.416, 0.936, 0.074, 0.938, 0.659, 0.035, 0.349, 0.169, 0.615, 0.464, 0.636, 0.173, 0.988, 0.172, 0.943, 0.223, 0.189, 0.021, 0.757, 0.126, 0.841, 0.154, 0.135, 0.37, 0.41, 0.685, 0.237, 0.828, 0.127, 0.543, 0.837, 0.119, 0.39, 0.719, 0.721, 0.793, 0.952, 0.7, 0.862, 0.887, 0.414, 0.161, 0.217, 0.183, 0.593, 0.581, 0.607, 0.923, 0.118, 0.804, 0.829, 0.68, 0.915, 0.008, 0.622, 0.926, 0.842, 0.583, 0.399, 0.063, 0.928, 0.437, 0.255, 0.233, 0.158, 0.775, 0.883, 0.923, 0.374, 0.45, 0.121, 0.754, 0.448, 0.465, 0.577, 0.2, 0.59, 0.038, 0.58, 0.492, 0.273, 0.145, 0.354, 0.978, 0.049, 0.508, 0.065, 0.845, 0.687, 0.57, 0.177, 0.175, 0.409, 0.588, 0.971, 0.608, 0.933, 0.74, 0.539, 0.51, 0.137]
global q = [0.875, 0.88, 0.979, 0.984, 0.661, 0.867, 0.506, 0.971, 0.715, 0.871, 0.976, 0.953, 0.558, 0.773, 0.944, 0.814, 0.413, 0.843, 0.316, 0.91, 0.981, 0.96, 0.403, 0.787, 0.964, 0.845, 0.987, 0.544, 0.783, 0.698, 0.781, 0.716, 0.855, 0.475, 0.957, 0.24, 0.873, 0.988, 0.959, 0.813, 0.765, 0.826, 0.752, 0.955, 0.237, 0.579, 0.986, 0.877, 0.848, 0.985, 0.993, 0.574, 0.833, 0.734, 0.745, 0.985, 0.367, 0.935, 0.774, 0.897, 0.767, 0.704, 0.969, 0.835, 0.338, 0.394, 0.85, 0.87, 0.61, 0.826, 0.984, 0.397, 0.883, 0.591, 0.838, 0.932, 0.811, 0.905, 0.958, 0.965, 0.689, 0.982, 0.743, 0.843, 0.565, 0.969, 0.783, 0.775, 0.992, 0.973, 0.763, 0.589, 0.616, 0.733, 0.417, 0.204, 0.974, 0.926, 0.834, 0.421, 0.798, 0.838, 0.429, 0.989, 0.7, 0.397, 0.97, 0.798, 0.982, 0.72, 0.857, 0.639, 0.93, 0.461, 0.813, 0.936, 0.894, 0.856, 0.803, 0.385, 0.273, 0.94, 0.955, 0.976, 0.334, 0.945, 0.867, 0.644, 0.997, 0.451, 0.557, 0.916, 0.987, 0.969, 0.26, 0.894, 0.898, 0.803, 0.538, 0.18, 0.422, 0.998, 0.984, 0.997, 0.926, 0.624, 0.464, 0.191, 0.962, 0.835, 0.895, 0.826, 0.264, 0.823, 0.507, 0.835, 0.59, 0.997, 0.779, 0.981, 0.878, 0.05, 0.57, 0.541, 0.957, 0.661, 0.856, 0.873, 0.999, 0.758, 0.974, 0.838, 0.522, 0.799, 0.807, 0.712, 0.979, 0.487, 0.82, 0.53, 0.782, 0.687, 0.426, 0.857, 0.76, 0.543, 0.876, 0.732, 0.486, 0.848, 0.897, 0.921, 0.987, 0.903, 0.942, 0.986, 0.75, 0.34, 0.863, 0.185, 0.991, 0.906, 0.798, 0.989, 0.553, 0.891, 0.89, 0.832, 0.955, 0.264, 0.761, 0.968, 0.963, 0.74, 0.951, 0.149, 0.994, 0.718, 0.282, 0.957, 0.375, 0.953, 0.939, 0.93, 0.773, 0.83, 0.276, 0.903, 0.755, 0.838, 0.909, 0.602, 0.782, 0.962, 0.781, 0.565, 0.617, 0.299, 0.355, 0.995, 0.251, 0.531, 0.098, 0.955, 0.737, 0.707, 0.402, 0.228, 0.534, 0.797, 0.974, 0.997, 0.996, 0.83, 0.87, 0.539, 0.853]
global origin = 1
global destination = 50