global arcs = [1 9; 1 10; 2 20; 2 34; 2 35; 3 10; 3 11; 3 24; 3 31; 4 2; 4 6; 4 23; 4 26; 4 28; 4 31; 5 3; 5 17; 5 26; 5 28; 6 3; 6 10; 6 12; 6 14; 6 27; 6 33; 7 31; 8 3; 8 35; 9 3; 9 13; 10 3; 10 24; 11 4; 11 5; 11 10; 11 25; 11 27; 11 31; 11 34; 12 9; 12 20; 12 22; 13 8; 13 10; 14 25; 14 31; 15 18; 15 30; 16 2; 16 4; 16 8; 16 15; 17 7; 17 23; 17 30; 17 31; 18 2; 18 12; 18 20; 18 29; 18 33; 19 22; 19 33; 19 35; 20 5; 20 11; 20 15; 20 24; 21 30; 22 10; 22 13; 22 31; 23 2; 23 12; 23 15; 23 24; 23 27; 24 13; 24 21; 24 26; 24 28; 24 32; 25 3; 25 32; 26 21; 26 33; 27 10; 28 15; 28 26; 28 27; 29 8; 29 10; 29 13; 30 7; 30 32; 31 2; 31 12; 31 29; 31 30; 32 5; 32 17; 32 34; 33 11; 33 16; 33 17; 33 19; 33 26; 34 8; 34 25]
global d_x = [4.0, 10.0, 8.0, 3.0, 3.0, 1.0, 5.0, 5.0, 9.0, 10.0, 7.0, 5.0, 6.0, 2.0, 2.0, 4.0, 5.0, 6.0, 2.0, 4.0, 4.0, 7.0, 7.0, 3.0, 4.0, 5.0, 1.0, 10.0, 8.0, 4.0, 9.0, 3.0, 3.0, 1.0, 2.0, 1.0, 6.0, 1.0, 1.0, 10.0, 3.0, 10.0, 8.0, 9.0, 9.0, 9.0, 9.0, 2.0, 2.0, 1.0, 4.0, 3.0, 9.0, 4.0, 9.0, 4.0, 10.0, 5.0, 7.0, 9.0, 6.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 3.0, 9.0, 7.0, 5.0, 10.0, 6.0, 8.0, 9.0, 6.0, 6.0, 7.0, 6.0, 3.0, 4.0, 9.0, 10.0, 4.0, 5.0, 8.0, 10.0, 8.0, 3.0, 4.0, 6.0, 9.0, 1.0, 3.0, 1.0, 1.0, 8.0, 7.0, 3.0, 1.0, 2.0, 9.0, 9.0, 9.0, 6.0, 7.0, 5.0, 8.0, 3.0]
global b_x = 5
global d_y = [2.0, 5.0, 10.0, 3.0, 3.0, 2.0, 7.0, 7.0, 2.0, 1.0, 7.0, 3.0, 7.0, 2.0, 7.0, 1.0, 1.0, 3.0, 3.0, 7.0, 1.0, 6.0, 5.0, 1.0, 10.0, 4.0, 2.0, 6.0, 4.0, 6.0, 7.0, 5.0, 9.0, 3.0, 10.0, 2.0, 1.0, 2.0, 6.0, 4.0, 5.0, 1.0, 3.0, 8.0, 1.0, 2.0, 8.0, 7.0, 3.0, 9.0, 8.0, 4.0, 4.0, 3.0, 8.0, 2.0, 6.0, 7.0, 5.0, 3.0, 10.0, 2.0, 4.0, 1.0, 1.0, 10.0, 9.0, 2.0, 3.0, 4.0, 4.0, 10.0, 5.0, 1.0, 7.0, 3.0, 6.0, 8.0, 9.0, 7.0, 2.0, 6.0, 4.0, 1.0, 8.0, 10.0, 10.0, 4.0, 4.0, 7.0, 3.0, 2.0, 3.0, 9.0, 7.0, 2.0, 2.0, 10.0, 7.0, 7.0, 4.0, 6.0, 9.0, 8.0, 8.0, 3.0, 7.0, 3.0, 10.0]
global b_y = 10
global p = [0.835, 0.459, 0.962, 0.112, 0.465, 0.701, 0.091, 0.687, 0.366, 0.664, 0.606, 0.976, 0.303, 0.863, 0.549, 0.882, 0.129, 0.406, 0.073, 0.144, 0.306, 0.817, 0.411, 0.627, 0.174, 0.676, 0.118, 0.66, 0.399, 0.453, 0.867, 0.57, 0.514, 0.429, 0.344, 0.741, 0.781, 0.803, 0.744, 0.583, 0.665, 0.366, 0.167, 0.612, 0.99, 0.329, 0.785, 0.738, 0.247, 0.223, 0.975, 0.318, 0.203, 0.267, 0.142, 0.96, 0.494, 0.078, 0.3, 0.907, 0.107, 0.469, 0.965, 0.305, 0.879, 0.38, 0.675, 0.638, 0.241, 0.104, 0.698, 0.285, 0.227, 0.093, 0.492, 0.307, 0.474, 0.105, 0.859, 0.83, 0.119, 0.44, 0.004, 0.754, 0.662, 0.458, 0.201, 0.884, 0.675, 0.497, 0.712, 0.31, 0.216, 0.403, 0.556, 0.902, 0.281, 0.862, 0.305, 0.958, 0.875, 0.661, 0.636, 0.605, 0.682, 0.557, 0.87, 0.073, 0.852]
global q = [0.859, 0.946, 0.999, 0.607, 0.654, 0.981, 0.208, 0.738, 0.539, 0.688, 0.957, 0.986, 0.452, 0.967, 0.94, 0.979, 0.388, 0.659, 0.991, 0.761, 0.535, 0.979, 0.647, 0.798, 0.683, 0.86, 0.404, 0.866, 0.425, 0.885, 0.888, 0.942, 0.67, 0.923, 0.404, 0.892, 0.923, 0.84, 0.799, 0.984, 0.986, 0.707, 0.867, 0.765, 0.993, 0.745, 0.982, 0.956, 0.826, 0.433, 0.991, 0.366, 0.219, 0.278, 0.147, 0.98, 0.763, 0.972, 0.397, 0.943, 0.277, 0.928, 0.968, 0.375, 0.993, 0.499, 0.691, 0.833, 0.318, 0.217, 0.994, 0.981, 0.929, 0.913, 0.695, 0.941, 0.5, 0.814, 0.866, 0.924, 0.557, 0.927, 0.643, 0.974, 0.98, 0.842, 0.918, 0.92, 0.997, 0.816, 0.832, 0.622, 0.275, 0.887, 0.65, 0.989, 0.5, 0.996, 0.902, 0.97, 0.962, 0.786, 0.998, 0.965, 0.732, 0.821, 0.904, 0.909, 0.915]
global origin = 1
global destination = 35