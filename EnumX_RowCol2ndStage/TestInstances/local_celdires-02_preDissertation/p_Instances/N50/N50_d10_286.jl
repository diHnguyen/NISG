global arcs = [1 2; 1 9; 1 19; 1 43; 1 50; 2 16; 2 23; 2 32; 2 35; 2 46; 2 49; 3 16; 3 22; 3 31; 3 40; 4 6; 4 9; 4 24; 4 26; 4 37; 4 41; 5 15; 5 29; 5 32; 5 44; 6 33; 6 34; 6 35; 7 36; 7 38; 7 42; 8 14; 8 17; 8 32; 8 35; 8 42; 9 3; 9 10; 9 24; 9 25; 9 47; 10 17; 10 34; 10 49; 10 50; 11 23; 11 27; 11 33; 11 36; 12 6; 12 25; 12 33; 12 49; 13 6; 13 12; 13 19; 13 41; 13 46; 13 48; 13 49; 13 50; 14 4; 14 10; 14 21; 14 26; 14 27; 14 47; 15 4; 15 17; 16 5; 16 6; 16 29; 16 33; 16 49; 17 4; 17 6; 17 12; 17 18; 17 24; 17 25; 17 26; 17 27; 17 39; 18 20; 18 33; 18 39; 18 49; 19 18; 19 26; 19 32; 19 36; 19 41; 19 42; 20 25; 20 47; 20 49; 20 50; 21 6; 21 12; 21 13; 21 30; 21 34; 21 45; 21 46; 22 19; 22 20; 22 26; 22 31; 22 47; 23 11; 23 30; 23 32; 23 33; 23 44; 24 6; 24 41; 25 16; 25 17; 25 18; 25 34; 26 10; 26 25; 26 47; 27 20; 27 41; 28 8; 28 13; 28 36; 28 37; 28 44; 28 46; 28 47; 29 5; 29 16; 29 23; 29 32; 29 36; 29 40; 29 41; 30 14; 30 39; 31 5; 31 25; 31 37; 31 49; 32 4; 32 21; 32 22; 32 30; 33 3; 33 22; 33 38; 33 45; 33 49; 34 10; 34 14; 34 28; 34 38; 35 5; 35 33; 35 40; 35 50; 36 10; 36 18; 36 21; 36 28; 36 31; 36 34; 36 47; 37 8; 37 16; 37 41; 37 42; 38 23; 38 26; 39 16; 39 43; 40 3; 40 18; 40 41; 41 8; 41 10; 41 28; 41 50; 42 18; 42 36; 42 46; 43 12; 43 24; 43 33; 43 46; 44 13; 44 32; 44 34; 44 37; 45 12; 45 14; 45 29; 46 2; 46 4; 46 10; 46 22; 46 26; 46 31; 46 47; 46 48; 47 3; 47 7; 47 16; 47 31; 47 36; 48 4; 48 13; 48 17; 48 21; 48 25; 48 27; 48 41; 48 46; 49 6; 49 12; 49 14; 49 32; 49 35; 49 40; 49 42]
global d_x = [9.0, 2.0, 7.0, 9.0, 10.0, 5.0, 6.0, 6.0, 2.0, 10.0, 4.0, 4.0, 9.0, 8.0, 9.0, 2.0, 5.0, 8.0, 4.0, 9.0, 8.0, 4.0, 7.0, 2.0, 1.0, 4.0, 1.0, 3.0, 3.0, 2.0, 4.0, 2.0, 8.0, 10.0, 1.0, 3.0, 9.0, 7.0, 6.0, 1.0, 3.0, 8.0, 6.0, 9.0, 8.0, 1.0, 2.0, 2.0, 4.0, 9.0, 6.0, 3.0, 10.0, 9.0, 1.0, 6.0, 3.0, 3.0, 6.0, 1.0, 5.0, 5.0, 9.0, 10.0, 8.0, 2.0, 10.0, 10.0, 1.0, 10.0, 8.0, 1.0, 1.0, 3.0, 10.0, 5.0, 10.0, 2.0, 2.0, 7.0, 6.0, 7.0, 9.0, 5.0, 4.0, 4.0, 8.0, 6.0, 10.0, 6.0, 3.0, 3.0, 8.0, 2.0, 1.0, 2.0, 3.0, 9.0, 6.0, 2.0, 1.0, 7.0, 10.0, 5.0, 6.0, 6.0, 8.0, 8.0, 2.0, 3.0, 10.0, 7.0, 10.0, 1.0, 6.0, 1.0, 4.0, 2.0, 3.0, 9.0, 3.0, 4.0, 8.0, 9.0, 9.0, 4.0, 10.0, 10.0, 6.0, 6.0, 10.0, 9.0, 4.0, 3.0, 1.0, 7.0, 8.0, 8.0, 10.0, 10.0, 7.0, 10.0, 3.0, 6.0, 2.0, 10.0, 8.0, 7.0, 8.0, 3.0, 10.0, 9.0, 2.0, 10.0, 6.0, 10.0, 7.0, 8.0, 4.0, 4.0, 2.0, 7.0, 3.0, 8.0, 6.0, 9.0, 4.0, 2.0, 9.0, 3.0, 7.0, 7.0, 5.0, 7.0, 10.0, 8.0, 8.0, 10.0, 7.0, 1.0, 3.0, 3.0, 3.0, 5.0, 4.0, 10.0, 9.0, 6.0, 6.0, 6.0, 1.0, 5.0, 1.0, 10.0, 7.0, 4.0, 4.0, 7.0, 6.0, 7.0, 9.0, 7.0, 9.0, 10.0, 10.0, 9.0, 9.0, 5.0, 10.0, 2.0, 5.0, 1.0, 7.0, 7.0, 5.0, 3.0, 9.0, 5.0, 1.0, 10.0, 3.0, 9.0, 5.0, 2.0, 6.0, 2.0]
global b_x = 5
global d_y = [3.0, 4.0, 10.0, 1.0, 6.0, 9.0, 2.0, 3.0, 3.0, 1.0, 7.0, 7.0, 6.0, 5.0, 6.0, 4.0, 6.0, 3.0, 5.0, 9.0, 1.0, 5.0, 6.0, 1.0, 1.0, 3.0, 5.0, 2.0, 9.0, 5.0, 10.0, 8.0, 10.0, 2.0, 10.0, 8.0, 10.0, 7.0, 6.0, 5.0, 1.0, 7.0, 10.0, 9.0, 3.0, 9.0, 9.0, 2.0, 3.0, 7.0, 1.0, 1.0, 3.0, 3.0, 9.0, 9.0, 10.0, 10.0, 5.0, 5.0, 5.0, 2.0, 7.0, 8.0, 1.0, 2.0, 6.0, 10.0, 4.0, 7.0, 4.0, 10.0, 4.0, 4.0, 9.0, 3.0, 3.0, 4.0, 4.0, 10.0, 8.0, 2.0, 2.0, 6.0, 9.0, 7.0, 7.0, 10.0, 4.0, 8.0, 6.0, 9.0, 6.0, 2.0, 4.0, 3.0, 6.0, 9.0, 2.0, 7.0, 10.0, 10.0, 9.0, 2.0, 10.0, 10.0, 10.0, 6.0, 2.0, 9.0, 10.0, 3.0, 6.0, 5.0, 7.0, 8.0, 4.0, 7.0, 10.0, 2.0, 9.0, 8.0, 10.0, 1.0, 4.0, 8.0, 1.0, 5.0, 5.0, 9.0, 4.0, 6.0, 4.0, 9.0, 5.0, 5.0, 2.0, 8.0, 9.0, 7.0, 9.0, 9.0, 4.0, 8.0, 1.0, 1.0, 1.0, 7.0, 4.0, 4.0, 2.0, 7.0, 5.0, 5.0, 2.0, 10.0, 7.0, 7.0, 9.0, 5.0, 2.0, 6.0, 8.0, 6.0, 6.0, 10.0, 1.0, 10.0, 1.0, 4.0, 1.0, 10.0, 7.0, 2.0, 7.0, 8.0, 2.0, 6.0, 5.0, 6.0, 10.0, 3.0, 7.0, 5.0, 8.0, 2.0, 5.0, 7.0, 7.0, 10.0, 9.0, 10.0, 3.0, 10.0, 8.0, 2.0, 9.0, 1.0, 6.0, 3.0, 8.0, 6.0, 10.0, 2.0, 10.0, 8.0, 2.0, 10.0, 8.0, 7.0, 6.0, 6.0, 4.0, 2.0, 10.0, 7.0, 5.0, 7.0, 6.0, 9.0, 3.0, 4.0, 7.0, 1.0, 8.0, 9.0]
global b_y = 10
global p = [0.53, 0.025, 0.51, 0.131, 0.045, 0.611, 0.645, 0.441, 0.469, 0.933, 0.171, 0.742, 0.99, 0.888, 0.893, 0.354, 0.376, 0.494, 0.414, 0.218, 0.838, 0.366, 0.585, 0.464, 0.607, 0.435, 0.804, 0.752, 0.916, 0.119, 0.187, 0.811, 0.479, 0.302, 0.029, 0.134, 0.927, 0.472, 0.045, 0.384, 0.486, 0.291, 0.675, 0.122, 0.108, 0.384, 0.357, 0.059, 0.034, 0.273, 0.736, 0.746, 0.551, 0.966, 0.007, 0.12, 0.456, 0.089, 0.476, 0.152, 0.293, 0.703, 0.263, 0.439, 0.173, 0.679, 0.318, 0.13, 0.14, 0.21, 0.745, 0.498, 0.796, 0.341, 0.583, 0.535, 0.586, 0.941, 0.937, 0.114, 0.773, 0.953, 0.556, 0.049, 0.684, 0.57, 0.387, 0.115, 0.213, 0.802, 0.011, 0.225, 0.636, 0.822, 0.446, 0.486, 0.756, 0.841, 0.124, 0.227, 0.961, 0.904, 0.732, 0.438, 0.958, 0.941, 0.206, 0.112, 0.952, 0.951, 0.254, 0.335, 0.308, 0.847, 0.775, 0.748, 0.311, 0.774, 0.46, 0.752, 0.869, 0.241, 0.977, 0.435, 0.098, 0.223, 0.428, 0.749, 0.582, 0.031, 0.519, 0.095, 0.98, 0.487, 0.1, 0.057, 0.186, 0.502, 0.097, 0.541, 0.862, 0.99, 0.678, 0.474, 0.442, 0.384, 0.756, 0.879, 0.294, 0.056, 0.407, 0.47, 0.955, 0.975, 0.002, 0.406, 0.829, 0.346, 0.198, 0.781, 0.087, 0.492, 0.034, 0.755, 0.472, 0.704, 0.103, 0.195, 0.624, 0.453, 0.026, 0.49, 0.579, 0.845, 0.11, 0.174, 0.614, 0.338, 0.019, 0.519, 0.046, 0.969, 0.673, 0.424, 0.133, 0.011, 0.478, 0.187, 0.173, 0.088, 0.334, 0.081, 0.532, 0.544, 0.577, 0.911, 0.592, 0.622, 0.725, 0.652, 0.513, 0.516, 0.389, 0.576, 0.528, 0.21, 0.141, 0.977, 0.004, 0.129, 0.818, 0.345, 0.456, 0.567, 0.062, 0.862, 0.544, 0.267, 0.27, 0.686, 0.356, 0.446, 0.308, 0.338, 0.75, 0.328]
global q = [0.557, 0.588, 0.926, 0.464, 0.72, 0.855, 0.961, 0.77, 0.73, 0.995, 0.748, 0.909, 0.991, 0.942, 0.96, 0.505, 0.49, 0.616, 0.678, 0.319, 0.875, 0.88, 0.679, 0.498, 0.997, 0.921, 0.936, 0.792, 0.972, 0.791, 0.262, 0.962, 0.619, 0.371, 0.091, 0.18, 0.933, 0.858, 0.428, 0.631, 0.632, 0.671, 0.699, 0.966, 0.348, 0.85, 0.499, 0.35, 0.889, 0.615, 0.891, 0.906, 0.863, 0.97, 0.569, 0.513, 0.533, 0.954, 0.508, 0.772, 0.301, 0.788, 0.484, 0.963, 0.746, 0.935, 0.901, 0.135, 0.545, 0.758, 0.877, 0.853, 0.845, 0.934, 0.901, 0.852, 0.881, 0.982, 0.939, 0.772, 0.917, 0.985, 0.692, 0.346, 0.757, 0.623, 0.85, 0.904, 0.556, 0.871, 0.836, 0.607, 0.942, 0.842, 0.682, 0.953, 0.879, 0.97, 0.517, 0.408, 0.983, 0.974, 0.943, 0.869, 0.978, 0.979, 0.946, 0.515, 0.958, 0.957, 0.812, 0.465, 0.84, 0.889, 0.945, 0.798, 0.4, 0.956, 0.92, 0.801, 0.999, 0.298, 0.983, 0.648, 0.53, 0.775, 0.516, 0.76, 0.877, 0.74, 0.84, 0.268, 0.997, 0.694, 0.719, 0.784, 0.933, 0.904, 0.324, 0.845, 0.981, 0.996, 0.73, 0.674, 0.964, 0.931, 0.977, 0.952, 0.975, 0.17, 0.93, 0.853, 0.994, 0.988, 0.013, 0.828, 0.9, 0.373, 0.591, 0.862, 0.58, 0.558, 0.989, 0.876, 0.788, 0.735, 0.271, 0.756, 0.72, 0.504, 0.136, 0.589, 0.653, 0.859, 0.334, 0.653, 0.696, 0.596, 0.532, 0.773, 0.31, 0.981, 0.799, 0.772, 0.135, 0.551, 0.719, 0.373, 0.475, 0.956, 0.515, 0.469, 0.84, 0.729, 0.594, 0.935, 0.99, 0.639, 0.817, 0.826, 0.715, 0.569, 0.793, 0.579, 0.948, 0.414, 0.789, 0.986, 0.233, 0.801, 0.894, 0.677, 0.856, 0.884, 0.364, 0.967, 0.925, 0.818, 0.711, 0.762, 0.64, 0.783, 0.399, 0.855, 0.904, 0.359]
global origin = 1
global destination = 50