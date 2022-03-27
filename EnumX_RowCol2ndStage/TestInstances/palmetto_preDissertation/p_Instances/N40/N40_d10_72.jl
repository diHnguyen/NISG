global arcs = [1 5; 1 38; 1 40; 2 13; 2 26; 3 17; 3 23; 3 27; 3 39; 4 5; 4 6; 4 11; 4 32; 4 34; 5 8; 6 3; 6 14; 6 30; 7 8; 7 16; 7 19; 7 28; 7 38; 8 4; 8 9; 8 31; 8 33; 8 34; 9 23; 9 24; 10 12; 10 30; 11 14; 11 20; 11 23; 11 30; 11 37; 12 10; 12 13; 12 21; 13 6; 13 11; 13 16; 13 21; 13 29; 14 10; 14 22; 14 33; 15 3; 15 7; 15 16; 15 22; 15 30; 15 36; 16 7; 16 10; 16 12; 16 23; 16 29; 16 34; 16 37; 17 10; 17 34; 18 21; 18 29; 18 37; 19 13; 19 25; 20 18; 20 23; 20 26; 20 29; 20 35; 20 37; 21 23; 22 5; 22 18; 22 25; 22 28; 23 2; 23 8; 23 27; 24 3; 24 10; 24 16; 24 33; 25 3; 25 4; 25 9; 25 12; 25 24; 25 30; 25 33; 26 9; 26 39; 27 7; 28 13; 28 18; 28 35; 29 16; 29 19; 30 3; 30 19; 30 24; 30 27; 31 12; 31 28; 31 29; 31 32; 32 26; 32 30; 33 5; 33 9; 33 10; 33 17; 33 21; 33 39; 34 15; 34 21; 34 24; 34 27; 35 11; 35 17; 35 19; 35 22; 35 23; 35 30; 35 37; 36 16; 36 18; 36 31; 36 33; 36 37; 36 38; 37 7; 37 17; 37 36; 37 40; 38 8; 38 11; 38 13; 38 18; 38 22; 38 35; 38 37; 38 39; 39 6; 39 7; 39 10; 39 16; 39 26; 39 27]
global d_x = [8.0, 5.0, 9.0, 2.0, 7.0, 5.0, 4.0, 1.0, 10.0, 6.0, 3.0, 8.0, 8.0, 1.0, 6.0, 8.0, 8.0, 7.0, 6.0, 6.0, 4.0, 7.0, 5.0, 10.0, 1.0, 8.0, 10.0, 5.0, 8.0, 2.0, 8.0, 10.0, 4.0, 5.0, 10.0, 10.0, 4.0, 2.0, 4.0, 2.0, 8.0, 2.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 8.0, 3.0, 9.0, 3.0, 8.0, 5.0, 6.0, 3.0, 5.0, 4.0, 8.0, 3.0, 9.0, 5.0, 8.0, 6.0, 7.0, 1.0, 10.0, 3.0, 7.0, 2.0, 1.0, 5.0, 8.0, 7.0, 2.0, 5.0, 6.0, 4.0, 1.0, 6.0, 7.0, 10.0, 2.0, 10.0, 10.0, 6.0, 3.0, 8.0, 5.0, 1.0, 3.0, 5.0, 4.0, 1.0, 5.0, 6.0, 1.0, 6.0, 3.0, 6.0, 8.0, 1.0, 10.0, 6.0, 1.0, 8.0, 8.0, 5.0, 6.0, 8.0, 1.0, 5.0, 8.0, 5.0, 3.0, 9.0, 3.0, 2.0, 6.0, 9.0, 5.0, 9.0, 9.0, 5.0, 6.0, 5.0, 9.0, 1.0, 4.0, 9.0, 4.0, 1.0, 3.0, 9.0, 7.0, 8.0, 8.0, 7.0, 8.0, 10.0, 9.0, 4.0, 3.0, 5.0, 9.0, 4.0, 1.0, 10.0, 10.0, 1.0, 4.0, 4.0]
global b_x = 5
global d_y = [4.0, 7.0, 10.0, 5.0, 2.0, 8.0, 7.0, 5.0, 1.0, 9.0, 6.0, 10.0, 7.0, 4.0, 9.0, 10.0, 9.0, 7.0, 10.0, 2.0, 6.0, 4.0, 6.0, 9.0, 10.0, 1.0, 10.0, 2.0, 7.0, 8.0, 7.0, 7.0, 7.0, 6.0, 7.0, 9.0, 7.0, 1.0, 2.0, 7.0, 10.0, 9.0, 9.0, 3.0, 4.0, 10.0, 4.0, 10.0, 1.0, 7.0, 2.0, 7.0, 2.0, 3.0, 3.0, 5.0, 6.0, 8.0, 5.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0, 9.0, 7.0, 5.0, 4.0, 1.0, 10.0, 10.0, 1.0, 8.0, 1.0, 3.0, 2.0, 7.0, 2.0, 2.0, 4.0, 4.0, 6.0, 5.0, 3.0, 4.0, 2.0, 5.0, 10.0, 9.0, 1.0, 2.0, 10.0, 8.0, 4.0, 5.0, 10.0, 7.0, 10.0, 2.0, 6.0, 8.0, 6.0, 2.0, 1.0, 4.0, 9.0, 4.0, 2.0, 10.0, 4.0, 3.0, 2.0, 2.0, 6.0, 8.0, 2.0, 7.0, 9.0, 2.0, 4.0, 7.0, 6.0, 9.0, 8.0, 6.0, 7.0, 1.0, 7.0, 5.0, 7.0, 7.0, 5.0, 2.0, 10.0, 8.0, 2.0, 8.0, 6.0, 2.0, 7.0, 9.0, 3.0, 5.0, 7.0, 2.0, 5.0, 3.0, 2.0, 9.0, 5.0, 2.0]
global b_y = 10
global p = [0.511, 0.137, 0.599, 0.173, 0.593, 0.817, 0.294, 0.934, 0.137, 0.908, 0.425, 0.118, 0.133, 0.708, 0.163, 0.44, 0.297, 0.621, 0.492, 0.643, 0.841, 0.443, 0.525, 0.903, 0.451, 0.943, 0.91, 0.267, 0.173, 0.621, 0.753, 0.706, 0.584, 0.639, 0.066, 0.992, 0.604, 0.551, 0.426, 0.989, 0.607, 0.295, 0.155, 0.473, 0.264, 0.349, 0.072, 0.936, 0.298, 0.581, 0.635, 0.993, 0.504, 0.634, 0.377, 0.685, 0.571, 0.872, 0.574, 0.308, 0.48, 0.964, 0.335, 0.222, 0.674, 0.354, 0.915, 0.9, 0.984, 0.071, 0.022, 0.892, 0.61, 0.227, 0.993, 0.492, 0.101, 0.849, 0.769, 0.247, 0.235, 0.641, 0.024, 0.348, 0.172, 0.457, 0.539, 0.32, 0.36, 0.182, 0.886, 0.502, 0.52, 0.741, 0.344, 0.283, 0.175, 0.066, 0.064, 0.592, 0.26, 0.036, 0.341, 0.339, 0.794, 0.491, 0.162, 0.151, 0.792, 0.018, 0.493, 0.259, 0.732, 0.11, 0.996, 0.33, 0.004, 0.042, 0.182, 0.516, 0.505, 0.37, 0.346, 0.548, 0.188, 0.086, 0.281, 0.695, 0.153, 0.758, 0.953, 0.275, 0.819, 0.157, 0.132, 0.185, 0.216, 0.315, 0.774, 0.865, 0.036, 0.886, 0.914, 0.625, 0.953, 0.668, 0.327, 0.744, 0.353, 0.22, 0.225, 0.896]
global q = [0.715, 0.281, 0.859, 0.624, 0.854, 0.915, 0.849, 0.975, 0.19, 0.96, 0.688, 0.716, 0.366, 0.737, 0.841, 0.853, 0.974, 0.654, 0.795, 0.939, 0.873, 0.847, 0.911, 0.908, 0.591, 0.972, 0.957, 0.469, 0.238, 0.776, 0.835, 0.924, 0.635, 0.922, 0.319, 0.996, 0.985, 0.831, 0.912, 0.997, 0.967, 0.399, 0.184, 0.779, 0.318, 0.624, 0.636, 0.941, 0.46, 0.997, 0.671, 0.994, 0.981, 0.833, 0.753, 0.881, 0.798, 0.881, 0.923, 0.312, 0.789, 0.996, 0.354, 0.402, 0.961, 0.568, 0.922, 0.991, 0.99, 0.201, 0.531, 0.969, 0.936, 0.912, 0.996, 0.597, 0.598, 0.936, 0.927, 0.655, 0.32, 0.905, 0.212, 0.595, 0.537, 0.587, 0.598, 0.952, 0.566, 0.712, 0.983, 0.712, 0.786, 0.997, 0.357, 0.866, 0.948, 0.943, 0.902, 0.639, 0.641, 0.758, 0.629, 0.576, 0.903, 0.935, 0.477, 0.551, 0.925, 0.885, 0.591, 0.935, 0.962, 0.767, 0.997, 0.439, 0.817, 0.651, 0.494, 0.968, 0.526, 0.913, 0.382, 0.737, 0.851, 0.938, 0.617, 0.717, 0.509, 0.875, 0.968, 0.31, 0.912, 0.248, 0.367, 0.574, 0.593, 0.617, 0.808, 0.902, 0.274, 0.912, 0.937, 0.784, 0.968, 0.759, 0.938, 0.778, 0.389, 0.328, 0.859, 0.956]
global origin = 1
global destination = 40