global arcs = [1 9; 1 19; 1 22; 1 33; 2 4; 2 13; 2 21; 2 30; 2 34; 2 37; 3 4; 3 15; 3 27; 4 23; 4 36; 5 2; 5 4; 5 6; 5 10; 5 14; 5 19; 5 28; 5 38; 5 40; 6 7; 6 12; 6 29; 7 3; 8 5; 8 12; 8 18; 9 29; 9 32; 9 34; 10 4; 10 7; 10 15; 10 30; 11 3; 11 7; 11 9; 11 15; 11 32; 12 22; 12 37; 13 16; 13 25; 14 19; 14 26; 14 27; 14 32; 15 17; 15 25; 15 26; 16 5; 16 6; 16 11; 16 13; 16 14; 16 21; 16 32; 16 37; 16 38; 17 5; 17 22; 17 23; 18 2; 18 3; 18 16; 18 23; 19 11; 19 12; 19 24; 19 39; 20 10; 20 14; 20 15; 20 25; 20 30; 20 33; 20 34; 20 40; 21 6; 21 28; 21 31; 22 31; 23 21; 23 24; 23 33; 23 35; 23 36; 23 37; 24 10; 24 33; 25 5; 25 6; 25 15; 25 23; 25 34; 26 19; 26 31; 26 32; 27 4; 27 20; 27 21; 27 32; 28 5; 28 7; 28 20; 29 2; 29 3; 29 9; 29 20; 29 30; 30 13; 30 24; 31 8; 31 14; 31 37; 31 39; 32 5; 32 26; 32 28; 32 34; 33 3; 33 10; 33 23; 33 34; 34 10; 34 11; 34 13; 35 13; 35 23; 35 25; 35 30; 35 32; 36 19; 36 20; 36 21; 36 25; 36 27; 36 32; 37 2; 37 3; 37 23; 37 34; 38 21; 38 28; 39 14; 39 16; 39 35; 39 37]
global d_x = [8.0, 8.0, 6.0, 9.0, 3.0, 1.0, 3.0, 4.0, 10.0, 9.0, 5.0, 7.0, 10.0, 3.0, 6.0, 3.0, 1.0, 7.0, 1.0, 9.0, 4.0, 4.0, 8.0, 7.0, 8.0, 3.0, 6.0, 9.0, 10.0, 1.0, 6.0, 2.0, 3.0, 7.0, 6.0, 4.0, 1.0, 4.0, 1.0, 4.0, 5.0, 4.0, 8.0, 5.0, 2.0, 8.0, 8.0, 8.0, 1.0, 2.0, 8.0, 2.0, 4.0, 1.0, 5.0, 9.0, 5.0, 4.0, 8.0, 7.0, 2.0, 3.0, 6.0, 3.0, 2.0, 2.0, 5.0, 4.0, 10.0, 7.0, 4.0, 8.0, 7.0, 5.0, 2.0, 4.0, 2.0, 7.0, 6.0, 5.0, 2.0, 3.0, 10.0, 8.0, 8.0, 5.0, 2.0, 3.0, 5.0, 3.0, 2.0, 10.0, 10.0, 10.0, 8.0, 4.0, 5.0, 5.0, 3.0, 2.0, 1.0, 7.0, 4.0, 8.0, 9.0, 2.0, 5.0, 7.0, 8.0, 5.0, 9.0, 2.0, 3.0, 6.0, 6.0, 8.0, 10.0, 10.0, 7.0, 1.0, 3.0, 9.0, 3.0, 1.0, 8.0, 5.0, 3.0, 7.0, 5.0, 2.0, 8.0, 3.0, 5.0, 3.0, 10.0, 7.0, 8.0, 7.0, 6.0, 9.0, 1.0, 5.0, 6.0, 8.0, 2.0, 3.0, 10.0, 9.0, 10.0, 8.0, 2.0, 10.0]
global b_x = 5
global d_y = [7.0, 3.0, 10.0, 1.0, 1.0, 5.0, 2.0, 6.0, 9.0, 7.0, 3.0, 7.0, 1.0, 6.0, 1.0, 2.0, 9.0, 3.0, 7.0, 10.0, 5.0, 6.0, 9.0, 4.0, 3.0, 2.0, 2.0, 3.0, 10.0, 8.0, 6.0, 8.0, 10.0, 10.0, 8.0, 9.0, 6.0, 3.0, 4.0, 9.0, 9.0, 10.0, 3.0, 7.0, 5.0, 7.0, 8.0, 4.0, 3.0, 7.0, 9.0, 3.0, 6.0, 2.0, 2.0, 7.0, 6.0, 1.0, 5.0, 9.0, 8.0, 9.0, 6.0, 2.0, 6.0, 1.0, 1.0, 10.0, 2.0, 2.0, 5.0, 10.0, 7.0, 8.0, 1.0, 6.0, 7.0, 4.0, 10.0, 9.0, 3.0, 8.0, 8.0, 7.0, 7.0, 5.0, 3.0, 10.0, 6.0, 5.0, 9.0, 6.0, 5.0, 8.0, 5.0, 5.0, 3.0, 7.0, 1.0, 7.0, 3.0, 3.0, 8.0, 9.0, 1.0, 4.0, 3.0, 2.0, 3.0, 6.0, 5.0, 2.0, 2.0, 6.0, 9.0, 5.0, 5.0, 9.0, 7.0, 9.0, 8.0, 2.0, 9.0, 7.0, 7.0, 6.0, 5.0, 6.0, 6.0, 1.0, 3.0, 6.0, 7.0, 9.0, 8.0, 6.0, 10.0, 4.0, 4.0, 4.0, 9.0, 8.0, 1.0, 6.0, 2.0, 7.0, 1.0, 6.0, 8.0, 3.0, 2.0, 8.0]
global b_y = 10
global p = [0.22, 0.212, 0.22, 0.225, 0.704, 0.156, 0.765, 0.487, 0.016, 0.586, 0.53, 0.167, 0.754, 0.43, 0.294, 0.457, 0.55, 0.806, 0.525, 0.633, 0.963, 0.197, 0.339, 0.37, 0.081, 0.031, 0.625, 0.975, 0.742, 0.819, 0.712, 0.219, 0.684, 0.474, 0.023, 0.532, 0.343, 0.441, 0.678, 0.942, 0.66, 0.492, 0.978, 0.543, 0.739, 0.61, 0.872, 0.594, 0.001, 0.499, 0.591, 0.05, 0.096, 0.957, 0.172, 0.549, 0.509, 0.114, 0.494, 0.147, 0.968, 0.547, 0.635, 0.78, 0.071, 0.64, 0.292, 0.639, 0.362, 0.627, 0.395, 0.42, 0.292, 0.89, 0.885, 0.93, 0.707, 0.009, 0.41, 0.533, 0.435, 0.806, 0.864, 0.982, 0.368, 0.555, 0.404, 0.332, 0.469, 0.936, 0.136, 0.603, 0.357, 0.382, 0.295, 0.89, 0.951, 0.442, 0.505, 0.657, 0.283, 0.958, 0.709, 0.028, 0.161, 0.832, 0.584, 0.01, 0.091, 0.323, 0.798, 0.235, 0.756, 0.495, 0.067, 0.466, 0.601, 0.259, 0.382, 0.822, 0.609, 0.443, 0.095, 0.032, 0.348, 0.846, 0.356, 0.282, 0.388, 0.954, 0.692, 0.466, 0.73, 0.227, 0.573, 0.123, 0.45, 0.304, 0.603, 0.333, 0.955, 0.754, 0.265, 0.577, 0.3, 0.124, 0.45, 0.082, 0.473, 0.766, 0.16, 0.324]
global q = [0.665, 0.951, 0.799, 0.337, 0.897, 0.484, 0.862, 0.806, 0.349, 0.842, 0.846, 0.783, 0.811, 0.862, 0.411, 0.896, 0.831, 0.919, 0.827, 0.661, 0.965, 0.536, 0.686, 0.627, 0.969, 0.586, 0.949, 0.977, 0.843, 0.999, 0.728, 0.267, 0.718, 0.61, 0.5, 0.786, 0.928, 0.715, 0.811, 0.946, 0.823, 0.768, 0.985, 0.585, 0.909, 0.765, 0.942, 0.845, 0.249, 0.7, 0.838, 0.251, 0.227, 0.975, 0.365, 0.795, 0.884, 0.346, 0.58, 0.638, 0.985, 0.9, 0.829, 0.996, 0.429, 0.805, 0.887, 0.783, 0.907, 0.713, 0.515, 0.675, 0.803, 0.99, 0.994, 0.987, 0.983, 0.625, 0.822, 0.96, 0.44, 0.885, 0.887, 0.983, 0.97, 0.588, 0.797, 0.627, 0.822, 0.964, 0.359, 0.783, 0.901, 0.591, 0.69, 0.91, 0.989, 0.538, 0.934, 0.965, 0.291, 0.969, 0.914, 0.391, 0.644, 0.973, 0.767, 0.737, 0.808, 0.88, 0.954, 0.662, 0.878, 0.944, 0.207, 0.756, 0.921, 0.687, 0.734, 0.902, 0.616, 0.715, 0.131, 0.963, 0.778, 0.902, 0.623, 0.828, 0.59, 0.984, 0.817, 0.788, 0.827, 0.518, 0.951, 0.232, 0.691, 0.449, 0.992, 0.487, 0.996, 0.952, 0.659, 0.678, 0.365, 0.166, 0.944, 0.403, 0.77, 0.988, 0.964, 0.646]
global origin = 1
global destination = 40