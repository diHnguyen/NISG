global arcs = [1 7; 1 21; 1 27; 2 9; 2 28; 3 5; 3 7; 3 21; 4 8; 4 13; 4 16; 4 19; 5 16; 6 2; 6 5; 6 12; 6 16; 6 28; 6 30; 7 23; 7 25; 7 28; 8 3; 8 40; 9 2; 9 4; 9 11; 9 16; 9 17; 9 18; 9 29; 9 35; 9 37; 10 18; 10 20; 10 21; 10 22; 10 24; 11 10; 11 16; 11 40; 12 2; 12 6; 12 19; 12 26; 13 17; 13 19; 13 22; 14 22; 14 26; 14 27; 15 7; 15 37; 15 39; 16 8; 16 12; 16 36; 16 39; 17 14; 17 18; 17 20; 17 28; 17 29; 17 35; 18 4; 18 7; 19 12; 19 17; 19 25; 19 26; 19 31; 19 38; 20 14; 20 18; 20 34; 21 4; 21 7; 21 10; 21 14; 21 15; 21 27; 21 30; 21 31; 22 3; 22 8; 22 12; 22 27; 22 29; 22 34; 22 36; 23 5; 23 9; 23 26; 23 30; 24 8; 24 10; 24 33; 24 38; 25 7; 25 9; 25 16; 25 26; 25 31; 26 3; 26 9; 26 18; 26 36; 26 37; 27 7; 27 11; 27 30; 27 35; 27 39; 28 27; 28 34; 28 37; 29 24; 29 25; 30 6; 30 22; 30 34; 31 6; 31 7; 31 34; 32 6; 32 29; 33 17; 33 39; 34 11; 34 16; 34 27; 34 30; 35 32; 36 14; 36 39; 36 40; 37 4; 37 10; 37 14; 37 17; 38 2; 38 10; 38 13; 38 21; 38 22; 38 23; 38 31; 39 9; 39 27]
global d_x = [3.0, 6.0, 1.0, 6.0, 9.0, 1.0, 10.0, 1.0, 8.0, 7.0, 9.0, 7.0, 1.0, 3.0, 4.0, 4.0, 7.0, 7.0, 4.0, 1.0, 6.0, 1.0, 9.0, 2.0, 8.0, 3.0, 4.0, 6.0, 7.0, 1.0, 6.0, 9.0, 4.0, 3.0, 4.0, 4.0, 1.0, 6.0, 8.0, 5.0, 4.0, 5.0, 8.0, 5.0, 2.0, 7.0, 5.0, 6.0, 9.0, 3.0, 4.0, 5.0, 4.0, 5.0, 5.0, 2.0, 3.0, 4.0, 3.0, 1.0, 3.0, 2.0, 5.0, 6.0, 7.0, 10.0, 10.0, 9.0, 10.0, 9.0, 4.0, 2.0, 2.0, 8.0, 6.0, 1.0, 8.0, 4.0, 10.0, 10.0, 8.0, 8.0, 8.0, 5.0, 8.0, 6.0, 9.0, 8.0, 5.0, 2.0, 3.0, 6.0, 7.0, 1.0, 4.0, 1.0, 9.0, 10.0, 9.0, 10.0, 8.0, 4.0, 9.0, 6.0, 9.0, 10.0, 7.0, 10.0, 7.0, 9.0, 8.0, 5.0, 10.0, 3.0, 5.0, 7.0, 5.0, 2.0, 5.0, 8.0, 1.0, 1.0, 6.0, 3.0, 4.0, 4.0, 6.0, 10.0, 7.0, 9.0, 9.0, 5.0, 8.0, 8.0, 1.0, 6.0, 8.0, 4.0, 6.0, 8.0, 4.0, 10.0, 7.0, 10.0, 7.0, 1.0, 10.0, 9.0, 9.0]
global b_x = 5
global d_y = [6.0, 1.0, 8.0, 4.0, 1.0, 6.0, 9.0, 2.0, 10.0, 8.0, 9.0, 9.0, 10.0, 5.0, 10.0, 10.0, 8.0, 7.0, 2.0, 3.0, 7.0, 1.0, 6.0, 8.0, 9.0, 1.0, 10.0, 3.0, 7.0, 1.0, 1.0, 8.0, 6.0, 1.0, 7.0, 2.0, 2.0, 4.0, 6.0, 10.0, 1.0, 1.0, 6.0, 8.0, 8.0, 6.0, 4.0, 5.0, 4.0, 7.0, 9.0, 1.0, 9.0, 8.0, 1.0, 6.0, 6.0, 5.0, 6.0, 8.0, 4.0, 8.0, 8.0, 2.0, 2.0, 9.0, 4.0, 9.0, 8.0, 3.0, 8.0, 10.0, 5.0, 1.0, 1.0, 9.0, 3.0, 5.0, 7.0, 10.0, 1.0, 9.0, 3.0, 7.0, 9.0, 2.0, 2.0, 8.0, 7.0, 4.0, 10.0, 10.0, 1.0, 8.0, 2.0, 5.0, 1.0, 6.0, 3.0, 3.0, 4.0, 1.0, 9.0, 9.0, 3.0, 1.0, 5.0, 7.0, 3.0, 7.0, 3.0, 6.0, 6.0, 6.0, 1.0, 4.0, 4.0, 1.0, 2.0, 8.0, 3.0, 5.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 4.0, 10.0, 1.0, 2.0, 6.0, 1.0, 4.0, 8.0, 9.0, 6.0, 7.0, 10.0, 7.0, 10.0, 3.0, 10.0, 9.0, 9.0, 4.0, 9.0, 2.0]
global b_y = 10
global p = [0.8, 0.157, 0.698, 0.364, 0.21, 0.648, 0.816, 0.579, 0.681, 0.809, 0.337, 0.048, 0.302, 0.721, 0.291, 0.98, 0.03, 0.027, 0.503, 0.139, 0.631, 0.74, 0.511, 0.534, 0.701, 0.566, 0.238, 0.347, 0.892, 0.524, 0.967, 0.232, 0.736, 0.547, 0.217, 0.344, 0.743, 0.713, 0.594, 0.556, 0.347, 0.706, 0.001, 0.771, 0.962, 0.224, 0.552, 0.18, 0.401, 0.403, 0.956, 0.984, 0.297, 0.101, 0.432, 0.131, 0.593, 0.077, 0.706, 0.679, 0.512, 0.098, 0.489, 0.671, 0.288, 0.441, 0.496, 0.596, 0.734, 0.554, 0.311, 0.057, 0.879, 0.616, 0.543, 0.05, 0.505, 0.078, 0.278, 0.708, 0.445, 0.663, 0.37, 0.533, 0.912, 0.548, 0.697, 0.487, 0.183, 0.324, 0.323, 0.506, 0.033, 0.263, 0.769, 0.335, 0.337, 0.486, 0.008, 0.075, 0.754, 0.039, 0.292, 0.825, 0.804, 0.624, 0.993, 0.71, 0.421, 0.226, 0.396, 0.068, 0.273, 0.02, 0.102, 0.592, 0.3, 0.768, 0.315, 0.6, 0.377, 0.699, 0.017, 0.326, 0.064, 0.361, 0.553, 0.9, 0.227, 0.418, 0.019, 0.796, 0.755, 0.591, 0.237, 0.135, 0.336, 0.535, 0.711, 0.884, 0.277, 0.749, 0.599, 0.143, 0.864, 0.651, 0.388, 0.787, 0.973]
global q = [0.909, 0.889, 0.834, 0.774, 0.435, 0.958, 0.91, 0.747, 0.942, 0.83, 0.64, 0.778, 0.385, 0.825, 0.409, 0.983, 0.113, 0.546, 0.893, 0.682, 0.965, 0.955, 0.587, 0.558, 0.813, 0.735, 0.69, 0.869, 0.935, 0.874, 0.987, 0.849, 0.857, 0.916, 0.218, 0.834, 0.896, 0.883, 0.815, 0.964, 0.719, 0.92, 0.82, 0.774, 0.963, 0.728, 0.802, 0.942, 0.877, 0.681, 0.968, 0.985, 0.409, 0.309, 0.665, 0.586, 0.695, 0.148, 0.906, 0.866, 0.735, 0.953, 0.758, 0.901, 0.389, 0.969, 0.535, 0.898, 0.779, 0.901, 0.869, 0.655, 0.99, 0.917, 0.557, 0.993, 0.821, 0.095, 0.772, 0.872, 0.657, 0.728, 0.766, 0.997, 0.996, 0.922, 0.913, 0.953, 0.836, 0.349, 0.406, 0.849, 0.266, 0.37, 0.8, 0.568, 0.695, 0.969, 0.245, 0.766, 0.78, 0.275, 0.442, 0.972, 0.893, 0.792, 0.995, 0.72, 0.88, 0.281, 0.642, 0.114, 0.636, 0.937, 0.498, 0.744, 0.592, 0.922, 0.578, 0.698, 0.848, 0.789, 0.543, 0.522, 0.542, 0.575, 0.798, 0.932, 0.931, 0.48, 0.489, 0.974, 0.771, 0.818, 0.596, 0.528, 0.399, 0.833, 0.732, 0.956, 0.848, 0.843, 0.795, 0.914, 0.989, 0.951, 0.966, 0.971, 0.974]
global origin = 1
global destination = 40