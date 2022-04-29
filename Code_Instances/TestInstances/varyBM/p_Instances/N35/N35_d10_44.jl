global arcs = [1 2; 1 16; 1 21; 1 31; 2 5; 2 17; 2 28; 2 30; 3 13; 3 23; 3 27; 3 29; 3 31; 3 35; 4 2; 4 18; 4 20; 4 31; 5 3; 5 11; 6 2; 7 16; 7 19; 7 24; 7 25; 7 32; 8 4; 8 7; 8 17; 8 19; 9 2; 9 27; 10 25; 10 29; 11 7; 11 22; 11 35; 12 9; 12 24; 13 14; 13 18; 13 19; 13 25; 13 28; 14 27; 15 5; 15 21; 15 28; 15 29; 15 33; 16 15; 16 18; 16 24; 16 31; 17 4; 17 5; 17 8; 17 24; 18 7; 18 8; 18 10; 19 4; 19 8; 19 13; 19 22; 19 26; 19 27; 20 13; 20 19; 20 23; 21 6; 21 10; 21 31; 22 3; 22 6; 22 9; 22 16; 22 30; 23 8; 23 16; 23 24; 23 33; 24 15; 24 16; 24 26; 24 32; 24 34; 25 7; 25 10; 25 15; 25 19; 26 4; 26 5; 26 7; 26 20; 26 33; 27 18; 28 21; 28 24; 28 25; 28 35; 29 13; 30 12; 30 17; 30 28; 30 34; 31 8; 31 17; 31 19; 31 21; 31 25; 32 2; 32 9; 32 13; 32 29; 32 35; 33 6; 33 27; 33 28; 34 20; 34 25]
global d_x = [4.0, 3.0, 10.0, 7.0, 8.0, 4.0, 4.0, 6.0, 9.0, 8.0, 2.0, 4.0, 3.0, 5.0, 2.0, 10.0, 7.0, 6.0, 1.0, 3.0, 6.0, 3.0, 3.0, 2.0, 6.0, 10.0, 9.0, 5.0, 10.0, 6.0, 9.0, 10.0, 6.0, 4.0, 6.0, 7.0, 9.0, 7.0, 2.0, 7.0, 7.0, 3.0, 2.0, 3.0, 8.0, 2.0, 8.0, 3.0, 6.0, 5.0, 5.0, 1.0, 9.0, 5.0, 8.0, 7.0, 3.0, 10.0, 3.0, 3.0, 10.0, 5.0, 4.0, 7.0, 4.0, 9.0, 8.0, 4.0, 3.0, 10.0, 10.0, 8.0, 2.0, 10.0, 1.0, 3.0, 2.0, 10.0, 10.0, 3.0, 4.0, 4.0, 6.0, 7.0, 4.0, 7.0, 3.0, 2.0, 8.0, 3.0, 1.0, 1.0, 10.0, 10.0, 6.0, 3.0, 4.0, 4.0, 8.0, 7.0, 9.0, 9.0, 6.0, 10.0, 6.0, 4.0, 4.0, 3.0, 1.0, 5.0, 8.0, 9.0, 8.0, 6.0, 5.0, 4.0, 1.0, 3.0, 7.0, 6.0, 6.0]
global b_x = 5
global d_y = [1.0, 7.0, 5.0, 5.0, 2.0, 1.0, 9.0, 4.0, 7.0, 9.0, 5.0, 4.0, 9.0, 9.0, 5.0, 8.0, 2.0, 7.0, 9.0, 3.0, 1.0, 7.0, 7.0, 9.0, 2.0, 8.0, 2.0, 3.0, 10.0, 1.0, 5.0, 1.0, 9.0, 6.0, 2.0, 4.0, 2.0, 3.0, 6.0, 3.0, 7.0, 9.0, 2.0, 9.0, 9.0, 4.0, 6.0, 7.0, 5.0, 2.0, 4.0, 6.0, 5.0, 6.0, 7.0, 9.0, 3.0, 7.0, 6.0, 3.0, 4.0, 7.0, 7.0, 9.0, 5.0, 8.0, 7.0, 4.0, 9.0, 6.0, 2.0, 1.0, 8.0, 7.0, 2.0, 3.0, 2.0, 5.0, 10.0, 5.0, 3.0, 6.0, 8.0, 7.0, 3.0, 2.0, 2.0, 4.0, 5.0, 5.0, 1.0, 5.0, 10.0, 2.0, 10.0, 7.0, 8.0, 7.0, 7.0, 8.0, 4.0, 2.0, 10.0, 9.0, 2.0, 5.0, 1.0, 3.0, 5.0, 2.0, 1.0, 5.0, 7.0, 5.0, 2.0, 3.0, 10.0, 1.0, 5.0, 5.0, 3.0]
global b_y = 10
global p = [0.075, 0.31, 0.596, 0.276, 0.692, 0.996, 0.426, 0.322, 0.659, 0.771, 0.744, 0.031, 0.399, 0.457, 0.403, 0.522, 0.265, 0.398, 0.702, 0.629, 0.96, 0.267, 0.854, 0.291, 0.474, 0.436, 0.297, 0.391, 0.046, 0.758, 0.729, 0.221, 0.352, 0.692, 0.22, 0.053, 0.642, 0.718, 0.394, 0.174, 0.921, 0.695, 0.346, 0.28, 0.061, 0.764, 0.046, 0.203, 0.226, 0.325, 0.981, 0.748, 0.897, 0.531, 0.285, 0.272, 0.237, 0.438, 0.728, 0.476, 0.851, 0.574, 0.175, 0.085, 0.934, 0.701, 0.414, 0.358, 0.366, 0.234, 0.967, 0.545, 0.372, 0.369, 0.999, 0.123, 0.724, 0.684, 0.438, 0.197, 0.034, 0.904, 0.791, 0.218, 0.072, 0.757, 0.95, 0.989, 0.096, 0.073, 0.033, 0.881, 0.43, 0.489, 0.58, 0.459, 0.073, 0.881, 0.623, 0.005, 0.584, 0.408, 0.599, 0.296, 0.211, 0.553, 0.97, 0.058, 0.761, 0.398, 0.664, 0.281, 0.296, 0.142, 0.306, 0.444, 0.475, 0.814, 0.935, 0.343, 0.332]
global q = [0.638, 0.922, 0.613, 0.383, 0.998, 0.996, 0.509, 0.743, 0.687, 0.778, 0.749, 0.69, 0.702, 0.738, 0.548, 0.987, 0.747, 0.509, 0.789, 0.689, 0.964, 0.799, 0.955, 0.796, 0.514, 0.62, 0.633, 0.762, 0.728, 0.997, 0.803, 0.343, 0.4, 0.735, 0.408, 0.2, 0.892, 0.818, 0.431, 0.318, 0.98, 0.968, 0.402, 0.427, 0.441, 0.858, 0.875, 0.915, 0.472, 0.328, 0.995, 0.831, 0.984, 0.581, 0.434, 0.911, 0.778, 0.96, 0.938, 0.653, 0.863, 0.846, 0.549, 0.9, 0.96, 0.882, 0.655, 0.472, 0.726, 0.35, 0.977, 0.859, 0.912, 0.54, 0.999, 0.266, 0.918, 0.963, 0.972, 0.271, 0.418, 0.914, 0.912, 0.724, 0.33, 0.796, 0.976, 0.997, 0.776, 0.328, 0.38, 0.917, 0.591, 0.892, 0.853, 0.798, 0.26, 0.884, 0.904, 0.178, 0.785, 0.547, 0.665, 0.755, 0.317, 0.678, 0.983, 0.869, 0.822, 0.963, 0.908, 0.888, 0.788, 0.35, 0.988, 0.766, 0.641, 0.924, 0.955, 0.527, 0.591]
global origin = 1
global destination = 35