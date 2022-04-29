global arcs = [1 10; 1 11; 1 13; 1 20; 1 24; 2 14; 2 35; 3 2; 4 3; 5 17; 5 33; 6 3; 6 19; 6 33; 7 2; 7 18; 8 10; 8 12; 8 18; 8 24; 9 2; 9 3; 9 10; 9 19; 9 22; 10 2; 10 8; 10 22; 11 18; 11 27; 12 19; 12 34; 13 6; 13 24; 13 30; 14 2; 14 11; 14 18; 14 20; 14 23; 14 26; 14 28; 15 17; 15 19; 16 13; 16 17; 16 26; 16 35; 17 27; 17 33; 18 5; 18 7; 18 9; 18 10; 18 21; 19 11; 19 14; 19 32; 19 34; 20 3; 20 8; 20 10; 20 16; 20 17; 20 19; 20 22; 21 11; 21 25; 21 26; 21 33; 22 20; 22 21; 22 27; 22 32; 23 3; 23 11; 24 3; 24 26; 24 27; 24 29; 25 15; 25 17; 25 21; 25 33; 26 18; 26 19; 27 3; 27 8; 27 9; 27 18; 27 21; 27 26; 28 9; 28 16; 28 21; 28 23; 28 31; 28 34; 29 4; 29 5; 29 13; 29 21; 29 22; 29 27; 29 28; 29 30; 29 34; 29 35; 30 2; 31 4; 31 8; 31 15; 31 25; 31 26; 32 10; 32 11; 33 6; 33 10; 34 2; 34 10; 34 13; 34 16; 34 24; 34 33]
global d_x = [8.0, 4.0, 10.0, 7.0, 6.0, 2.0, 10.0, 6.0, 3.0, 1.0, 9.0, 5.0, 6.0, 6.0, 8.0, 5.0, 7.0, 9.0, 7.0, 3.0, 8.0, 9.0, 1.0, 10.0, 6.0, 2.0, 3.0, 1.0, 10.0, 4.0, 5.0, 5.0, 5.0, 6.0, 9.0, 10.0, 10.0, 8.0, 8.0, 1.0, 6.0, 7.0, 10.0, 8.0, 9.0, 9.0, 9.0, 2.0, 1.0, 6.0, 9.0, 8.0, 6.0, 2.0, 1.0, 2.0, 2.0, 2.0, 8.0, 4.0, 6.0, 4.0, 2.0, 9.0, 1.0, 2.0, 10.0, 5.0, 6.0, 3.0, 7.0, 5.0, 3.0, 7.0, 7.0, 2.0, 9.0, 7.0, 8.0, 9.0, 3.0, 3.0, 1.0, 8.0, 9.0, 2.0, 1.0, 1.0, 7.0, 10.0, 6.0, 8.0, 2.0, 2.0, 7.0, 3.0, 7.0, 9.0, 6.0, 3.0, 4.0, 5.0, 3.0, 3.0, 5.0, 2.0, 10.0, 9.0, 1.0, 10.0, 6.0, 8.0, 1.0, 10.0, 9.0, 1.0, 5.0, 3.0, 10.0, 6.0, 8.0, 9.0, 8.0, 10.0]
global b_x = 5
global d_y = [6.0, 10.0, 9.0, 1.0, 4.0, 9.0, 4.0, 5.0, 6.0, 4.0, 6.0, 2.0, 1.0, 1.0, 10.0, 4.0, 9.0, 10.0, 8.0, 2.0, 3.0, 2.0, 2.0, 3.0, 9.0, 3.0, 6.0, 7.0, 3.0, 3.0, 8.0, 7.0, 3.0, 4.0, 7.0, 5.0, 7.0, 6.0, 9.0, 4.0, 8.0, 10.0, 4.0, 8.0, 8.0, 4.0, 10.0, 9.0, 4.0, 3.0, 1.0, 1.0, 10.0, 9.0, 5.0, 5.0, 8.0, 9.0, 5.0, 7.0, 5.0, 8.0, 1.0, 10.0, 8.0, 10.0, 7.0, 8.0, 2.0, 10.0, 4.0, 3.0, 8.0, 8.0, 3.0, 2.0, 9.0, 2.0, 10.0, 10.0, 4.0, 7.0, 10.0, 10.0, 4.0, 1.0, 9.0, 9.0, 6.0, 9.0, 8.0, 6.0, 5.0, 3.0, 8.0, 4.0, 8.0, 6.0, 2.0, 3.0, 6.0, 4.0, 7.0, 6.0, 1.0, 10.0, 6.0, 7.0, 6.0, 6.0, 10.0, 1.0, 7.0, 5.0, 4.0, 5.0, 10.0, 4.0, 10.0, 4.0, 8.0, 9.0, 3.0, 6.0]
global b_y = 10
global p = [0.288, 0.865, 0.309, 0.807, 0.758, 0.266, 0.064, 0.998, 0.782, 0.594, 0.081, 0.063, 0.502, 0.774, 0.905, 0.694, 0.527, 0.09, 0.002, 0.749, 0.489, 0.104, 0.572, 0.22, 0.601, 0.813, 0.231, 0.48, 0.98, 0.819, 0.376, 0.843, 0.117, 0.1, 0.864, 0.898, 0.66, 0.009, 0.704, 0.697, 0.422, 0.22, 0.486, 0.665, 0.079, 0.811, 0.62, 0.931, 0.424, 0.693, 0.857, 0.438, 0.81, 0.008, 0.479, 0.351, 0.996, 0.727, 0.619, 0.556, 0.31, 0.731, 0.563, 0.874, 0.39, 0.012, 0.867, 0.396, 0.618, 0.076, 0.86, 0.269, 0.081, 0.286, 0.187, 0.734, 0.377, 0.061, 0.332, 0.357, 0.667, 0.757, 0.112, 0.937, 0.033, 0.783, 0.342, 0.783, 0.242, 0.963, 0.634, 0.279, 0.773, 0.092, 0.604, 0.917, 0.706, 0.622, 0.763, 0.955, 0.337, 0.549, 0.61, 0.966, 0.59, 0.351, 0.51, 0.355, 0.113, 0.518, 0.393, 0.852, 0.832, 0.29, 0.224, 0.091, 0.399, 0.794, 0.095, 0.746, 0.867, 0.454, 0.265, 0.632]
global q = [0.632, 0.951, 0.373, 0.998, 0.937, 0.37, 0.273, 0.998, 0.9, 0.685, 0.434, 0.763, 0.788, 0.862, 0.907, 0.916, 0.701, 0.435, 0.282, 0.997, 0.931, 0.285, 0.638, 0.866, 0.615, 0.906, 0.661, 0.891, 0.985, 0.907, 0.446, 0.925, 0.345, 0.103, 0.92, 0.964, 0.678, 0.984, 0.871, 0.731, 0.703, 0.233, 0.811, 0.794, 0.665, 0.97, 0.726, 0.944, 0.899, 0.725, 0.919, 0.83, 0.884, 0.576, 0.487, 0.693, 0.998, 0.796, 0.91, 0.755, 0.564, 0.995, 0.883, 0.99, 0.855, 0.63, 0.984, 0.42, 0.73, 0.365, 0.979, 0.283, 0.95, 0.738, 0.515, 0.891, 0.985, 0.543, 0.743, 0.704, 0.758, 0.789, 0.256, 0.99, 0.885, 0.862, 0.697, 0.895, 0.359, 0.986, 0.929, 0.496, 0.898, 0.345, 0.901, 0.977, 0.86, 0.674, 0.875, 0.99, 0.558, 0.839, 0.631, 0.993, 0.726, 0.839, 0.789, 0.573, 0.199, 0.77, 0.896, 0.932, 0.885, 0.314, 0.562, 0.916, 0.683, 0.879, 0.493, 0.908, 0.986, 0.73, 0.618, 0.65]
global origin = 1
global destination = 35