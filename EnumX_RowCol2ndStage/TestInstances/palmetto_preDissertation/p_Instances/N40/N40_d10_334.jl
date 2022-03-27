global arcs = [1 2; 1 19; 1 23; 1 26; 1 37; 1 39; 1 40; 2 6; 2 10; 2 14; 2 15; 2 16; 2 26; 2 35; 2 39; 3 7; 3 28; 3 40; 4 2; 4 29; 5 6; 5 12; 5 20; 5 38; 6 3; 6 11; 6 14; 6 29; 7 3; 7 6; 7 23; 7 28; 7 33; 8 11; 8 18; 8 26; 8 28; 8 29; 8 31; 8 38; 9 16; 9 26; 9 36; 10 6; 10 13; 10 31; 11 7; 11 10; 11 20; 11 31; 11 33; 12 11; 13 8; 13 15; 13 18; 13 20; 13 32; 14 6; 14 24; 15 2; 15 4; 15 9; 16 8; 16 21; 16 26; 16 36; 17 8; 17 10; 17 26; 17 39; 18 2; 18 10; 18 13; 18 23; 18 32; 19 2; 19 13; 19 31; 19 40; 20 2; 20 16; 20 37; 20 38; 21 2; 21 8; 21 12; 21 17; 21 19; 21 27; 21 31; 22 13; 22 15; 22 40; 23 19; 23 21; 23 39; 24 3; 24 7; 24 34; 24 35; 24 40; 25 16; 25 27; 26 15; 26 38; 26 39; 27 11; 27 13; 27 18; 28 19; 28 20; 28 39; 29 5; 29 8; 29 13; 29 23; 29 28; 29 38; 29 39; 30 2; 30 5; 30 9; 30 12; 30 31; 31 7; 31 15; 31 22; 31 28; 31 29; 32 11; 32 14; 33 10; 33 19; 33 22; 34 18; 34 22; 34 28; 34 31; 34 32; 34 33; 34 35; 34 37; 35 7; 35 21; 35 27; 36 22; 36 25; 36 33; 36 34; 37 2; 37 36; 38 2; 38 7; 38 10; 38 24; 38 31; 38 37; 39 30; 39 32]
global d_x = [2.0, 9.0, 8.0, 7.0, 9.0, 2.0, 7.0, 3.0, 1.0, 10.0, 1.0, 3.0, 10.0, 5.0, 7.0, 2.0, 2.0, 3.0, 1.0, 1.0, 3.0, 8.0, 4.0, 8.0, 4.0, 3.0, 8.0, 5.0, 4.0, 1.0, 4.0, 7.0, 7.0, 7.0, 7.0, 1.0, 9.0, 9.0, 4.0, 2.0, 4.0, 4.0, 6.0, 2.0, 5.0, 6.0, 10.0, 2.0, 8.0, 8.0, 1.0, 6.0, 10.0, 6.0, 3.0, 7.0, 4.0, 5.0, 8.0, 6.0, 5.0, 7.0, 5.0, 7.0, 2.0, 1.0, 9.0, 2.0, 1.0, 1.0, 5.0, 6.0, 7.0, 1.0, 5.0, 7.0, 7.0, 7.0, 9.0, 10.0, 6.0, 6.0, 2.0, 2.0, 4.0, 2.0, 9.0, 4.0, 9.0, 5.0, 2.0, 6.0, 7.0, 7.0, 3.0, 6.0, 6.0, 4.0, 8.0, 3.0, 7.0, 2.0, 9.0, 1.0, 4.0, 9.0, 7.0, 7.0, 7.0, 10.0, 3.0, 3.0, 6.0, 8.0, 3.0, 7.0, 4.0, 9.0, 2.0, 10.0, 1.0, 7.0, 4.0, 1.0, 7.0, 2.0, 8.0, 7.0, 10.0, 5.0, 5.0, 7.0, 1.0, 4.0, 3.0, 3.0, 4.0, 2.0, 10.0, 6.0, 6.0, 9.0, 3.0, 9.0, 5.0, 5.0, 6.0, 7.0, 9.0, 3.0, 9.0, 9.0, 1.0, 4.0, 5.0, 10.0, 3.0, 4.0, 8.0]
global b_x = 5
global d_y = [10.0, 5.0, 1.0, 6.0, 7.0, 7.0, 7.0, 2.0, 5.0, 1.0, 1.0, 7.0, 1.0, 2.0, 6.0, 1.0, 6.0, 2.0, 3.0, 1.0, 7.0, 3.0, 8.0, 9.0, 1.0, 7.0, 8.0, 6.0, 6.0, 10.0, 7.0, 2.0, 9.0, 7.0, 5.0, 7.0, 5.0, 7.0, 3.0, 10.0, 4.0, 7.0, 3.0, 8.0, 5.0, 5.0, 4.0, 9.0, 9.0, 2.0, 4.0, 9.0, 2.0, 6.0, 4.0, 7.0, 8.0, 7.0, 8.0, 4.0, 5.0, 2.0, 9.0, 9.0, 10.0, 1.0, 8.0, 4.0, 7.0, 7.0, 4.0, 8.0, 10.0, 2.0, 3.0, 7.0, 7.0, 6.0, 6.0, 5.0, 9.0, 6.0, 9.0, 1.0, 4.0, 6.0, 8.0, 4.0, 3.0, 10.0, 1.0, 7.0, 2.0, 4.0, 8.0, 2.0, 4.0, 10.0, 4.0, 2.0, 4.0, 1.0, 10.0, 8.0, 2.0, 9.0, 5.0, 7.0, 1.0, 2.0, 5.0, 2.0, 5.0, 4.0, 4.0, 5.0, 8.0, 9.0, 6.0, 1.0, 2.0, 7.0, 3.0, 4.0, 7.0, 2.0, 2.0, 9.0, 5.0, 4.0, 4.0, 3.0, 2.0, 9.0, 8.0, 3.0, 4.0, 2.0, 9.0, 5.0, 6.0, 6.0, 10.0, 1.0, 6.0, 4.0, 8.0, 1.0, 4.0, 5.0, 9.0, 7.0, 2.0, 4.0, 4.0, 8.0, 9.0, 3.0, 7.0]
global b_y = 10
global p = [0.22, 0.774, 0.398, 0.607, 0.193, 0.078, 0.258, 0.087, 0.092, 0.467, 0.863, 0.971, 0.493, 0.289, 0.704, 0.914, 0.224, 0.293, 0.621, 0.967, 0.636, 0.482, 0.259, 0.093, 0.923, 0.161, 0.926, 0.215, 0.415, 0.397, 0.193, 0.423, 0.678, 0.24, 0.308, 0.046, 0.458, 0.954, 0.888, 0.544, 0.798, 0.344, 0.018, 0.557, 0.646, 0.371, 0.399, 0.258, 0.652, 0.127, 0.014, 0.899, 0.867, 0.2, 0.915, 0.954, 0.555, 0.97, 0.481, 0.061, 0.702, 0.202, 0.388, 0.523, 0.881, 0.639, 0.807, 0.826, 0.351, 0.471, 0.351, 0.274, 0.82, 0.137, 0.049, 0.229, 0.007, 0.686, 0.259, 0.416, 0.032, 0.242, 0.52, 0.152, 0.133, 0.778, 0.667, 0.448, 0.756, 0.428, 0.095, 0.09, 0.092, 0.367, 0.221, 0.197, 0.052, 0.669, 0.971, 0.377, 0.39, 0.419, 0.017, 0.998, 0.271, 0.949, 0.021, 0.955, 0.785, 0.227, 0.671, 0.667, 0.67, 0.524, 0.729, 0.428, 0.856, 0.901, 0.791, 0.762, 0.056, 0.224, 0.488, 0.724, 0.175, 0.205, 0.818, 0.553, 0.988, 0.865, 0.592, 0.165, 0.242, 0.752, 0.427, 0.457, 0.386, 0.921, 0.23, 0.814, 0.571, 0.251, 0.693, 0.283, 0.227, 0.449, 0.64, 0.718, 0.63, 0.028, 0.382, 0.702, 0.589, 0.675, 0.447, 0.391, 0.83, 0.421, 0.281]
global q = [0.577, 0.994, 0.426, 0.948, 0.512, 0.87, 0.88, 0.329, 0.231, 0.925, 0.874, 0.989, 0.544, 0.976, 0.743, 0.972, 0.555, 0.953, 0.64, 0.972, 0.93, 0.867, 0.764, 0.206, 0.947, 0.548, 0.969, 0.946, 0.557, 0.905, 0.971, 0.676, 0.928, 0.255, 0.85, 0.28, 0.607, 0.958, 0.932, 0.568, 0.993, 0.376, 0.987, 0.72, 0.681, 0.863, 0.542, 0.832, 0.762, 0.816, 0.542, 0.954, 0.924, 0.298, 0.954, 0.981, 0.642, 0.998, 0.984, 0.674, 0.735, 0.574, 0.874, 0.891, 0.925, 0.907, 0.911, 0.998, 0.727, 0.774, 0.575, 0.9, 0.907, 0.539, 0.253, 0.413, 0.325, 0.713, 0.824, 0.635, 0.974, 0.869, 0.774, 0.586, 0.341, 0.906, 0.925, 0.943, 0.997, 0.98, 0.126, 0.426, 0.675, 0.729, 0.792, 0.362, 0.056, 0.695, 0.99, 0.735, 0.672, 0.957, 0.096, 0.998, 0.391, 0.994, 0.228, 0.985, 0.98, 0.714, 0.83, 0.92, 0.883, 0.923, 0.878, 0.82, 0.926, 0.931, 0.96, 0.894, 0.683, 0.368, 0.509, 0.956, 0.67, 0.505, 0.894, 0.656, 0.995, 0.867, 0.653, 0.566, 0.476, 0.925, 0.614, 0.476, 0.715, 0.981, 0.424, 0.873, 0.799, 0.937, 0.982, 0.788, 0.861, 0.452, 0.812, 0.729, 0.65, 0.358, 0.993, 0.888, 0.986, 0.799, 0.466, 0.908, 0.854, 0.741, 0.705]
global origin = 1
global destination = 40