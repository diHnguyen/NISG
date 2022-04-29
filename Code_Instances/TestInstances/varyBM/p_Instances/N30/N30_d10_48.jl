global arcs = [1 6; 1 21; 2 5; 2 10; 2 19; 2 20; 2 23; 3 2; 3 14; 3 15; 3 22; 4 25; 5 10; 5 11; 5 17; 6 5; 6 13; 6 18; 6 23; 7 15; 7 30; 8 15; 8 22; 9 4; 10 5; 10 24; 10 29; 10 30; 11 2; 11 16; 11 27; 12 9; 12 22; 13 6; 13 11; 13 25; 14 16; 14 17; 15 8; 15 9; 16 21; 17 2; 17 3; 17 8; 17 11; 17 12; 17 18; 17 26; 18 2; 18 9; 18 24; 19 2; 19 16; 20 2; 20 9; 20 11; 20 25; 20 28; 21 3; 21 6; 21 7; 21 11; 21 24; 22 11; 22 12; 22 21; 22 28; 23 6; 23 14; 23 20; 23 22; 24 17; 24 22; 25 11; 25 29; 26 2; 26 11; 26 12; 26 16; 26 21; 27 16; 28 4; 28 11; 28 14; 28 15; 28 17; 28 20; 28 26; 29 4; 29 13; 29 23]
global d_x = [3.0, 9.0, 1.0, 3.0, 7.0, 5.0, 1.0, 5.0, 2.0, 6.0, 5.0, 5.0, 1.0, 7.0, 5.0, 2.0, 9.0, 1.0, 5.0, 10.0, 4.0, 2.0, 2.0, 3.0, 4.0, 7.0, 8.0, 8.0, 9.0, 5.0, 9.0, 6.0, 8.0, 2.0, 7.0, 8.0, 6.0, 3.0, 4.0, 5.0, 9.0, 6.0, 2.0, 6.0, 4.0, 7.0, 3.0, 8.0, 4.0, 5.0, 4.0, 2.0, 9.0, 9.0, 6.0, 8.0, 5.0, 8.0, 5.0, 8.0, 3.0, 4.0, 8.0, 2.0, 5.0, 6.0, 10.0, 1.0, 4.0, 9.0, 2.0, 9.0, 2.0, 10.0, 10.0, 1.0, 5.0, 6.0, 4.0, 9.0, 5.0, 3.0, 3.0, 5.0, 8.0, 2.0, 7.0, 8.0, 3.0, 4.0, 9.0]
global b_x = 5
global d_y = [5.0, 3.0, 3.0, 4.0, 3.0, 1.0, 9.0, 2.0, 3.0, 2.0, 1.0, 2.0, 3.0, 3.0, 10.0, 8.0, 4.0, 3.0, 10.0, 7.0, 10.0, 1.0, 8.0, 9.0, 3.0, 1.0, 10.0, 9.0, 5.0, 6.0, 2.0, 8.0, 6.0, 5.0, 9.0, 9.0, 6.0, 6.0, 9.0, 2.0, 2.0, 9.0, 10.0, 10.0, 1.0, 2.0, 7.0, 1.0, 10.0, 10.0, 7.0, 10.0, 1.0, 10.0, 3.0, 6.0, 6.0, 10.0, 9.0, 3.0, 8.0, 9.0, 10.0, 5.0, 6.0, 1.0, 5.0, 8.0, 7.0, 6.0, 7.0, 10.0, 4.0, 6.0, 1.0, 1.0, 8.0, 1.0, 9.0, 2.0, 5.0, 9.0, 2.0, 3.0, 8.0, 5.0, 5.0, 6.0, 8.0, 1.0, 4.0]
global b_y = 10
global p = [0.956, 0.767, 0.33, 0.297, 0.655, 0.804, 0.105, 0.046, 0.411, 0.139, 0.72, 0.797, 0.996, 0.036, 0.443, 0.059, 0.843, 0.155, 0.736, 0.645, 0.418, 0.557, 0.955, 0.846, 0.896, 0.811, 0.966, 0.241, 0.188, 0.321, 0.815, 0.204, 0.987, 0.093, 0.774, 0.787, 0.084, 0.834, 0.151, 0.988, 0.924, 0.114, 0.252, 0.506, 0.621, 0.236, 0.514, 0.352, 0.571, 0.127, 0.208, 0.562, 0.99, 0.339, 0.052, 0.79, 0.919, 0.44, 0.872, 0.943, 0.163, 0.219, 0.667, 0.414, 0.965, 0.525, 0.048, 0.414, 0.212, 0.711, 0.398, 0.579, 0.766, 0.985, 0.858, 0.761, 0.018, 0.541, 0.682, 0.705, 0.536, 0.26, 0.09, 0.943, 0.813, 0.799, 0.074, 0.435, 0.737, 0.682, 0.654]
global q = [0.983, 0.801, 0.938, 0.97, 0.702, 0.932, 0.384, 0.111, 0.884, 0.794, 0.773, 0.882, 0.998, 0.312, 0.642, 0.972, 0.983, 0.993, 0.9, 0.743, 0.508, 0.738, 0.956, 0.922, 0.9, 0.908, 0.985, 0.659, 0.685, 0.85, 0.875, 0.465, 0.992, 0.412, 0.904, 0.917, 0.864, 0.986, 0.883, 0.995, 0.975, 0.311, 0.401, 0.779, 0.727, 0.876, 0.751, 0.603, 0.901, 0.326, 0.623, 0.802, 0.996, 0.791, 0.222, 0.8, 0.988, 0.667, 0.962, 0.962, 0.825, 0.987, 0.982, 0.724, 0.993, 0.73, 0.818, 0.891, 0.585, 0.826, 0.414, 0.793, 0.778, 0.99, 0.941, 0.98, 0.483, 0.934, 0.765, 0.895, 0.663, 0.344, 0.537, 0.961, 0.948, 0.939, 0.374, 0.626, 0.885, 0.828, 0.656]
global origin = 1
global destination = 30