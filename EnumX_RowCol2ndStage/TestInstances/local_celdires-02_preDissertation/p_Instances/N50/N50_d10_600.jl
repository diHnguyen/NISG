global arcs = [1 3; 1 8; 1 12; 1 13; 1 18; 1 25; 1 27; 1 32; 1 38; 1 46; 1 49; 2 10; 2 36; 2 46; 3 7; 3 27; 3 30; 4 5; 4 8; 4 21; 5 11; 5 16; 5 28; 5 33; 5 34; 5 39; 5 44; 5 50; 6 9; 6 43; 6 48; 6 50; 7 14; 7 23; 7 41; 7 43; 8 2; 8 6; 8 11; 8 22; 8 24; 8 43; 9 5; 9 17; 9 44; 9 47; 10 2; 10 14; 10 20; 11 2; 11 41; 12 10; 12 36; 12 38; 13 12; 13 16; 13 24; 13 49; 14 16; 14 32; 15 3; 15 7; 15 13; 15 38; 15 42; 16 5; 16 11; 16 30; 16 37; 16 39; 16 48; 17 32; 17 35; 17 50; 18 11; 18 42; 19 5; 19 9; 19 46; 19 49; 19 50; 20 2; 20 5; 20 17; 20 23; 20 37; 21 32; 21 38; 21 39; 21 41; 22 17; 22 20; 22 24; 22 28; 22 32; 22 45; 23 10; 23 20; 23 39; 23 43; 24 15; 24 33; 24 34; 25 12; 25 17; 25 19; 25 20; 25 26; 25 43; 25 44; 25 46; 26 12; 26 15; 26 23; 26 32; 26 35; 26 38; 26 49; 27 20; 27 41; 27 50; 28 3; 28 6; 28 37; 29 33; 30 5; 30 24; 30 28; 30 35; 30 38; 30 46; 31 17; 31 29; 31 47; 32 13; 32 17; 32 37; 32 43; 32 49; 33 2; 33 18; 33 34; 33 40; 33 48; 34 3; 34 16; 34 21; 34 31; 34 39; 35 10; 35 16; 35 23; 35 37; 35 50; 36 23; 36 27; 37 4; 37 11; 37 29; 37 42; 38 2; 38 4; 38 10; 38 14; 38 33; 38 37; 38 45; 39 22; 39 30; 39 43; 40 11; 40 13; 40 21; 40 22; 40 25; 40 33; 41 8; 41 9; 41 27; 42 5; 42 8; 42 15; 42 19; 42 26; 42 29; 42 32; 42 36; 42 40; 42 45; 42 50; 43 4; 43 6; 43 18; 43 19; 43 40; 44 12; 44 16; 44 25; 44 39; 44 47; 45 7; 45 9; 45 11; 45 14; 45 23; 45 31; 45 42; 46 6; 46 33; 47 10; 47 39; 47 41; 47 42; 48 10; 48 16; 48 20; 48 24; 48 27; 48 44; 48 45; 48 50; 49 4; 49 7; 49 9; 49 22; 49 33; 49 42]
global d_x = [6.0, 5.0, 2.0, 3.0, 3.0, 6.0, 2.0, 7.0, 6.0, 9.0, 4.0, 10.0, 5.0, 8.0, 10.0, 6.0, 5.0, 5.0, 3.0, 6.0, 5.0, 5.0, 5.0, 2.0, 2.0, 3.0, 3.0, 8.0, 4.0, 5.0, 7.0, 10.0, 8.0, 1.0, 7.0, 10.0, 10.0, 3.0, 2.0, 6.0, 1.0, 10.0, 2.0, 5.0, 6.0, 4.0, 10.0, 3.0, 5.0, 5.0, 9.0, 9.0, 10.0, 8.0, 2.0, 1.0, 1.0, 10.0, 4.0, 8.0, 9.0, 5.0, 9.0, 9.0, 8.0, 8.0, 2.0, 3.0, 1.0, 3.0, 9.0, 1.0, 1.0, 6.0, 6.0, 5.0, 8.0, 4.0, 7.0, 10.0, 6.0, 10.0, 3.0, 5.0, 4.0, 8.0, 10.0, 2.0, 1.0, 1.0, 10.0, 8.0, 2.0, 8.0, 4.0, 2.0, 2.0, 1.0, 6.0, 7.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 7.0, 1.0, 9.0, 3.0, 8.0, 5.0, 9.0, 4.0, 5.0, 5.0, 4.0, 10.0, 2.0, 1.0, 3.0, 8.0, 8.0, 7.0, 8.0, 5.0, 7.0, 7.0, 10.0, 3.0, 10.0, 1.0, 1.0, 7.0, 3.0, 2.0, 10.0, 9.0, 7.0, 3.0, 10.0, 5.0, 5.0, 3.0, 1.0, 10.0, 7.0, 2.0, 8.0, 8.0, 9.0, 2.0, 6.0, 9.0, 6.0, 9.0, 2.0, 6.0, 9.0, 5.0, 8.0, 2.0, 5.0, 5.0, 7.0, 6.0, 4.0, 8.0, 8.0, 9.0, 5.0, 2.0, 10.0, 2.0, 5.0, 8.0, 3.0, 4.0, 1.0, 1.0, 3.0, 2.0, 8.0, 4.0, 5.0, 10.0, 3.0, 3.0, 6.0, 8.0, 7.0, 5.0, 6.0, 4.0, 4.0, 5.0, 4.0, 9.0, 7.0, 3.0, 3.0, 6.0, 6.0, 9.0, 4.0, 5.0, 10.0, 7.0, 8.0, 3.0, 1.0, 4.0, 4.0, 7.0, 8.0, 7.0, 9.0, 6.0, 10.0, 9.0, 5.0, 9.0, 1.0, 3.0, 5.0, 6.0, 1.0]
global b_x = 5
global d_y = [9.0, 6.0, 1.0, 6.0, 2.0, 6.0, 5.0, 2.0, 9.0, 1.0, 1.0, 5.0, 10.0, 5.0, 1.0, 1.0, 2.0, 9.0, 4.0, 9.0, 10.0, 3.0, 5.0, 10.0, 10.0, 10.0, 1.0, 9.0, 8.0, 7.0, 9.0, 6.0, 10.0, 7.0, 10.0, 4.0, 4.0, 3.0, 3.0, 4.0, 1.0, 2.0, 8.0, 2.0, 8.0, 2.0, 5.0, 3.0, 1.0, 5.0, 9.0, 2.0, 5.0, 8.0, 5.0, 7.0, 5.0, 10.0, 3.0, 9.0, 2.0, 7.0, 8.0, 3.0, 8.0, 9.0, 9.0, 8.0, 9.0, 7.0, 9.0, 1.0, 10.0, 3.0, 5.0, 6.0, 2.0, 1.0, 5.0, 2.0, 4.0, 9.0, 9.0, 4.0, 9.0, 8.0, 10.0, 3.0, 2.0, 7.0, 7.0, 1.0, 2.0, 4.0, 1.0, 3.0, 3.0, 6.0, 10.0, 5.0, 4.0, 4.0, 1.0, 7.0, 2.0, 7.0, 7.0, 1.0, 9.0, 6.0, 5.0, 9.0, 6.0, 6.0, 4.0, 10.0, 2.0, 9.0, 3.0, 4.0, 10.0, 10.0, 3.0, 10.0, 6.0, 4.0, 10.0, 10.0, 10.0, 8.0, 9.0, 5.0, 1.0, 4.0, 2.0, 9.0, 10.0, 5.0, 8.0, 3.0, 8.0, 1.0, 4.0, 7.0, 3.0, 10.0, 8.0, 10.0, 9.0, 8.0, 3.0, 10.0, 1.0, 8.0, 3.0, 2.0, 9.0, 1.0, 5.0, 4.0, 7.0, 9.0, 10.0, 4.0, 6.0, 6.0, 9.0, 8.0, 7.0, 2.0, 4.0, 3.0, 9.0, 4.0, 7.0, 1.0, 3.0, 5.0, 8.0, 2.0, 5.0, 2.0, 7.0, 5.0, 5.0, 1.0, 7.0, 5.0, 1.0, 9.0, 10.0, 8.0, 1.0, 7.0, 10.0, 5.0, 3.0, 6.0, 6.0, 7.0, 1.0, 4.0, 4.0, 1.0, 6.0, 8.0, 3.0, 8.0, 6.0, 4.0, 5.0, 9.0, 3.0, 6.0, 1.0, 3.0, 1.0, 9.0, 3.0, 8.0, 6.0, 9.0, 7.0, 5.0, 4.0, 3.0, 9.0]
global b_y = 10
global p = [0.139, 0.575, 0.216, 0.424, 0.888, 0.623, 0.857, 0.535, 0.391, 0.272, 0.408, 0.655, 0.981, 0.841, 0.832, 0.137, 0.069, 0.639, 0.584, 0.915, 0.604, 0.1, 0.392, 0.566, 0.217, 0.683, 0.535, 0.974, 0.808, 0.269, 0.296, 0.504, 0.502, 0.838, 0.323, 0.096, 0.806, 0.4, 0.418, 0.383, 0.978, 0.047, 0.666, 0.092, 0.348, 0.273, 0.593, 0.947, 0.746, 0.109, 0.744, 0.843, 0.172, 0.932, 0.805, 0.626, 0.792, 0.797, 0.693, 0.302, 0.579, 0.024, 0.288, 0.158, 0.372, 0.996, 0.343, 0.273, 0.04, 0.645, 0.628, 0.529, 0.735, 0.989, 0.721, 0.492, 0.297, 0.21, 0.041, 0.533, 0.189, 0.251, 0.884, 0.908, 0.621, 0.052, 0.679, 0.418, 0.502, 0.565, 0.286, 0.162, 0.146, 0.129, 0.57, 0.905, 0.066, 0.072, 0.143, 0.451, 0.692, 0.508, 0.726, 0.644, 0.889, 0.527, 0.202, 0.07, 0.645, 0.707, 0.055, 0.097, 0.454, 0.767, 0.697, 0.47, 0.183, 0.231, 0.823, 0.803, 0.292, 0.224, 0.564, 0.401, 0.63, 0.145, 0.333, 0.262, 0.925, 0.485, 0.087, 0.77, 0.149, 0.011, 0.821, 0.993, 0.941, 0.971, 0.086, 0.444, 0.186, 0.707, 0.442, 0.447, 0.221, 0.677, 0.541, 0.369, 0.551, 0.08, 0.218, 0.015, 0.722, 0.852, 0.519, 0.371, 0.539, 0.365, 0.32, 0.05, 0.478, 0.471, 0.781, 0.407, 0.927, 0.574, 0.095, 0.176, 0.971, 0.757, 0.227, 0.508, 0.592, 0.725, 0.366, 0.174, 0.411, 0.943, 0.801, 0.69, 0.635, 0.258, 0.512, 0.897, 0.609, 0.191, 0.238, 0.172, 0.188, 0.784, 0.01, 0.191, 0.685, 0.433, 0.404, 0.663, 0.204, 0.381, 0.525, 0.121, 0.68, 0.209, 0.66, 0.553, 0.856, 0.038, 0.613, 0.58, 0.903, 0.681, 0.28, 0.061, 0.282, 0.701, 0.275, 0.579, 0.639, 0.621, 0.315, 0.834, 0.254, 0.068, 0.838, 0.824, 0.326, 0.043, 0.742]
global q = [0.691, 0.794, 0.587, 0.932, 0.933, 0.636, 0.999, 0.917, 0.746, 0.293, 0.955, 0.763, 0.992, 0.936, 0.946, 0.878, 0.357, 0.792, 0.951, 0.964, 0.758, 0.671, 0.471, 0.894, 0.908, 0.742, 0.607, 0.996, 0.897, 0.798, 0.695, 0.831, 0.989, 0.866, 0.8, 0.928, 0.879, 0.608, 0.685, 0.512, 0.994, 0.973, 0.949, 0.884, 0.436, 0.729, 0.593, 0.992, 0.795, 0.594, 0.952, 0.999, 0.503, 0.976, 0.855, 0.957, 0.812, 0.816, 0.937, 0.488, 0.823, 0.627, 0.538, 0.342, 0.462, 0.997, 0.382, 0.674, 0.084, 0.647, 0.723, 0.619, 0.82, 0.996, 0.815, 0.786, 0.432, 0.29, 0.58, 0.742, 0.486, 0.714, 0.987, 0.993, 0.981, 0.318, 0.834, 0.451, 0.527, 0.732, 0.684, 0.772, 0.296, 0.808, 0.864, 0.931, 0.888, 0.703, 0.765, 0.804, 0.753, 0.995, 0.985, 0.918, 0.91, 0.564, 0.875, 0.442, 0.876, 0.829, 0.719, 0.454, 0.798, 0.91, 0.885, 0.962, 0.562, 0.669, 0.944, 0.807, 0.777, 0.983, 0.811, 0.503, 0.771, 0.359, 0.917, 0.629, 0.991, 0.768, 0.841, 0.892, 0.563, 0.634, 0.845, 0.993, 0.962, 0.987, 0.334, 0.872, 0.913, 0.978, 0.643, 0.729, 0.771, 0.999, 0.573, 0.391, 0.742, 0.242, 0.755, 0.901, 0.857, 0.894, 0.534, 0.481, 0.625, 0.399, 0.941, 0.558, 0.698, 0.524, 0.996, 0.879, 0.98, 0.69, 0.262, 0.74, 0.976, 0.789, 0.931, 0.583, 0.828, 0.777, 0.829, 0.3, 0.625, 0.962, 0.846, 0.779, 0.762, 0.587, 0.772, 0.963, 0.965, 0.676, 0.54, 0.714, 0.402, 0.843, 0.999, 0.531, 0.83, 0.672, 0.619, 0.858, 0.87, 0.955, 0.933, 0.97, 0.827, 0.741, 0.699, 0.625, 0.933, 0.289, 0.872, 0.68, 0.908, 0.78, 0.861, 0.097, 0.324, 0.976, 0.557, 0.991, 0.741, 0.964, 0.587, 0.903, 0.962, 0.912, 0.935, 0.997, 0.641, 0.938, 0.98]
global origin = 1
global destination = 50