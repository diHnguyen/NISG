global arcs = [1 5; 1 29; 1 41; 2 4; 2 6; 2 10; 2 12; 2 14; 2 16; 2 22; 2 24; 2 32; 2 35; 2 36; 3 4; 3 10; 3 14; 3 19; 3 35; 3 36; 3 38; 3 44; 4 34; 4 48; 5 4; 5 8; 5 17; 5 22; 5 26; 5 40; 5 41; 6 13; 6 20; 6 26; 6 31; 6 33; 6 36; 6 42; 6 44; 6 49; 7 14; 7 16; 7 18; 7 20; 7 28; 7 36; 7 39; 7 49; 8 3; 8 4; 8 14; 8 16; 8 30; 8 46; 9 11; 9 28; 9 29; 9 30; 9 32; 10 9; 10 23; 10 24; 10 25; 11 10; 11 14; 11 31; 11 35; 11 43; 11 45; 11 48; 12 8; 12 20; 12 25; 13 16; 13 17; 13 20; 13 24; 13 44; 13 47; 14 13; 14 29; 14 37; 14 46; 15 25; 15 43; 16 11; 16 19; 16 26; 16 33; 16 42; 16 49; 16 50; 17 4; 17 11; 17 21; 17 24; 17 40; 18 3; 18 22; 18 24; 18 43; 18 44; 19 18; 19 31; 19 32; 19 43; 20 4; 20 8; 20 25; 20 26; 20 31; 20 32; 21 13; 21 17; 21 19; 22 12; 22 16; 22 50; 23 10; 23 16; 23 42; 24 4; 24 10; 24 23; 24 27; 24 30; 24 41; 24 49; 25 7; 25 8; 25 18; 25 22; 25 27; 25 30; 25 35; 25 44; 25 45; 25 48; 26 10; 26 12; 26 14; 26 19; 26 29; 26 33; 26 35; 26 41; 26 49; 27 2; 27 3; 27 4; 27 8; 27 18; 27 30; 27 33; 27 36; 27 46; 28 6; 28 10; 28 16; 28 21; 28 30; 28 41; 28 47; 29 14; 29 21; 30 18; 30 39; 31 9; 31 26; 31 32; 32 19; 32 43; 33 7; 33 11; 33 13; 33 27; 33 30; 33 34; 34 42; 34 49; 34 50; 35 5; 35 7; 35 13; 35 16; 35 18; 35 27; 35 28; 36 11; 36 15; 36 28; 36 29; 36 33; 36 37; 36 38; 36 43; 36 49; 37 4; 37 5; 37 10; 37 11; 37 24; 38 9; 38 43; 38 45; 39 2; 39 20; 40 23; 40 27; 40 35; 40 39; 41 6; 41 11; 41 22; 41 23; 41 25; 41 36; 42 20; 42 27; 42 45; 43 7; 43 18; 43 22; 43 40; 43 41; 43 45; 44 3; 44 9; 44 22; 44 33; 44 38; 44 46; 44 47; 45 3; 45 5; 45 11; 45 15; 45 23; 45 29; 45 35; 46 14; 46 17; 46 31; 46 38; 47 48; 48 13; 48 16; 48 25; 48 29; 48 30; 48 36; 49 11; 49 36; 49 42; 49 43; 49 50]
global d_x = [8.0, 3.0, 5.0, 5.0, 1.0, 3.0, 7.0, 1.0, 5.0, 2.0, 4.0, 2.0, 7.0, 6.0, 6.0, 9.0, 7.0, 3.0, 6.0, 10.0, 10.0, 3.0, 8.0, 4.0, 10.0, 6.0, 6.0, 9.0, 2.0, 7.0, 3.0, 4.0, 7.0, 8.0, 3.0, 5.0, 7.0, 1.0, 6.0, 4.0, 3.0, 8.0, 9.0, 10.0, 7.0, 6.0, 4.0, 6.0, 5.0, 9.0, 6.0, 8.0, 2.0, 9.0, 6.0, 6.0, 10.0, 6.0, 1.0, 8.0, 10.0, 4.0, 4.0, 3.0, 1.0, 6.0, 8.0, 5.0, 5.0, 5.0, 5.0, 3.0, 9.0, 1.0, 8.0, 4.0, 6.0, 7.0, 5.0, 4.0, 7.0, 1.0, 10.0, 6.0, 5.0, 6.0, 8.0, 2.0, 7.0, 2.0, 2.0, 1.0, 6.0, 3.0, 8.0, 3.0, 7.0, 9.0, 4.0, 5.0, 6.0, 9.0, 7.0, 6.0, 1.0, 3.0, 6.0, 7.0, 2.0, 6.0, 6.0, 7.0, 7.0, 5.0, 6.0, 1.0, 1.0, 4.0, 1.0, 7.0, 6.0, 1.0, 7.0, 1.0, 3.0, 7.0, 10.0, 6.0, 1.0, 1.0, 3.0, 1.0, 3.0, 2.0, 6.0, 7.0, 9.0, 6.0, 9.0, 4.0, 4.0, 8.0, 5.0, 5.0, 6.0, 9.0, 4.0, 7.0, 9.0, 2.0, 5.0, 3.0, 4.0, 6.0, 7.0, 2.0, 10.0, 6.0, 8.0, 9.0, 1.0, 1.0, 8.0, 3.0, 7.0, 2.0, 8.0, 8.0, 2.0, 2.0, 3.0, 1.0, 8.0, 8.0, 1.0, 2.0, 2.0, 1.0, 2.0, 2.0, 9.0, 4.0, 10.0, 7.0, 3.0, 4.0, 10.0, 2.0, 5.0, 1.0, 6.0, 8.0, 2.0, 5.0, 3.0, 2.0, 9.0, 3.0, 8.0, 3.0, 7.0, 10.0, 6.0, 4.0, 3.0, 5.0, 1.0, 1.0, 9.0, 2.0, 7.0, 1.0, 3.0, 4.0, 8.0, 4.0, 7.0, 8.0, 1.0, 9.0, 9.0, 6.0, 8.0, 1.0, 3.0, 3.0, 10.0, 7.0, 7.0, 7.0, 7.0, 9.0, 2.0, 1.0, 8.0, 3.0, 4.0, 7.0, 8.0, 7.0, 4.0, 5.0, 10.0, 4.0, 9.0, 2.0, 2.0, 8.0, 2.0, 2.0, 9.0, 10.0, 4.0, 9.0, 6.0, 1.0]
global b_x = 5
global d_y = [5.0, 8.0, 3.0, 7.0, 1.0, 6.0, 4.0, 2.0, 1.0, 1.0, 6.0, 6.0, 4.0, 6.0, 4.0, 8.0, 6.0, 8.0, 10.0, 3.0, 9.0, 6.0, 9.0, 10.0, 4.0, 8.0, 6.0, 2.0, 2.0, 6.0, 6.0, 6.0, 5.0, 3.0, 10.0, 9.0, 7.0, 5.0, 7.0, 4.0, 1.0, 7.0, 5.0, 6.0, 4.0, 2.0, 10.0, 10.0, 10.0, 9.0, 3.0, 10.0, 3.0, 5.0, 10.0, 5.0, 8.0, 6.0, 10.0, 4.0, 4.0, 7.0, 6.0, 10.0, 6.0, 5.0, 1.0, 5.0, 5.0, 2.0, 10.0, 2.0, 2.0, 10.0, 6.0, 4.0, 8.0, 2.0, 7.0, 7.0, 9.0, 10.0, 1.0, 10.0, 10.0, 1.0, 1.0, 3.0, 3.0, 9.0, 2.0, 2.0, 4.0, 5.0, 2.0, 10.0, 9.0, 3.0, 4.0, 7.0, 8.0, 8.0, 2.0, 2.0, 8.0, 2.0, 2.0, 6.0, 5.0, 8.0, 5.0, 8.0, 1.0, 10.0, 2.0, 1.0, 3.0, 6.0, 8.0, 7.0, 8.0, 2.0, 1.0, 4.0, 7.0, 6.0, 9.0, 7.0, 8.0, 10.0, 5.0, 5.0, 9.0, 5.0, 8.0, 8.0, 2.0, 3.0, 10.0, 9.0, 10.0, 8.0, 6.0, 6.0, 10.0, 2.0, 6.0, 7.0, 8.0, 9.0, 10.0, 6.0, 2.0, 10.0, 4.0, 4.0, 1.0, 4.0, 6.0, 4.0, 9.0, 7.0, 3.0, 10.0, 1.0, 3.0, 2.0, 3.0, 3.0, 6.0, 6.0, 3.0, 2.0, 7.0, 2.0, 2.0, 7.0, 6.0, 2.0, 3.0, 5.0, 10.0, 8.0, 6.0, 2.0, 2.0, 3.0, 5.0, 7.0, 4.0, 9.0, 4.0, 6.0, 5.0, 3.0, 6.0, 2.0, 1.0, 6.0, 7.0, 5.0, 2.0, 8.0, 9.0, 8.0, 5.0, 2.0, 1.0, 8.0, 9.0, 3.0, 3.0, 9.0, 1.0, 7.0, 1.0, 8.0, 6.0, 3.0, 10.0, 2.0, 2.0, 9.0, 6.0, 1.0, 7.0, 5.0, 8.0, 7.0, 5.0, 1.0, 10.0, 4.0, 7.0, 4.0, 9.0, 7.0, 9.0, 1.0, 2.0, 10.0, 8.0, 10.0, 9.0, 9.0, 2.0, 7.0, 4.0, 3.0, 4.0, 4.0, 6.0, 3.0, 2.0, 3.0, 8.0]
global b_y = 10
global p = [0.435, 0.226, 0.36, 0.704, 0.156, 0.375, 0.139, 0.177, 0.507, 0.032, 0.296, 0.789, 0.811, 0.345, 0.448, 0.867, 0.497, 0.102, 0.552, 0.571, 0.242, 0.721, 0.355, 0.774, 0.021, 0.652, 0.973, 0.646, 0.767, 0.867, 0.159, 0.126, 0.163, 0.943, 0.944, 0.175, 0.281, 0.572, 0.982, 0.362, 0.561, 0.385, 0.244, 0.727, 0.307, 0.481, 0.362, 0.249, 0.742, 0.134, 0.015, 0.919, 0.651, 0.02, 0.24, 0.807, 0.285, 0.269, 0.721, 0.334, 0.044, 0.896, 0.81, 0.589, 0.998, 0.376, 0.071, 0.007, 0.86, 0.178, 0.063, 0.981, 0.995, 0.586, 0.57, 0.513, 0.548, 0.309, 0.166, 0.272, 0.592, 0.997, 0.335, 0.983, 0.981, 0.628, 0.304, 0.091, 0.261, 0.227, 0.704, 0.615, 0.333, 0.825, 0.6, 0.914, 0.97, 0.951, 0.417, 0.253, 0.928, 0.597, 0.127, 0.834, 0.22, 0.205, 0.311, 0.508, 0.235, 0.954, 0.381, 0.893, 0.375, 0.862, 0.398, 0.571, 0.073, 0.689, 0.382, 0.991, 0.078, 0.47, 0.021, 0.992, 0.007, 0.325, 0.835, 0.203, 0.713, 0.163, 0.462, 0.688, 0.818, 0.971, 0.475, 0.04, 0.04, 0.752, 0.017, 0.832, 0.808, 0.194, 0.041, 0.69, 0.651, 0.644, 0.2, 0.982, 0.177, 0.938, 0.392, 0.384, 0.175, 0.472, 0.077, 0.473, 0.11, 0.975, 0.065, 0.964, 0.311, 0.86, 0.893, 0.791, 0.19, 0.432, 0.587, 0.618, 0.991, 0.132, 0.411, 0.925, 0.956, 0.865, 0.665, 0.224, 0.396, 0.581, 0.399, 0.575, 0.592, 0.704, 0.308, 0.778, 0.762, 0.021, 0.747, 0.8, 0.376, 0.454, 0.081, 0.536, 0.948, 0.677, 0.045, 0.429, 0.191, 0.469, 0.497, 0.655, 0.039, 0.081, 0.706, 0.546, 0.006, 0.979, 0.546, 0.754, 0.762, 0.025, 0.853, 0.793, 0.672, 0.313, 0.714, 0.484, 0.465, 0.287, 0.221, 0.782, 0.441, 0.833, 0.372, 0.978, 0.884, 0.012, 0.322, 0.534, 0.876, 0.069, 0.439, 0.162, 0.423, 0.077, 0.992, 0.651, 0.338, 0.616, 0.329, 0.477, 0.498, 0.071, 0.156, 0.133, 0.223, 0.928, 0.081, 0.645, 0.031, 0.672, 0.378, 0.432, 0.818, 0.905, 0.712, 0.852]
global q = [0.976, 0.502, 0.963, 0.83, 0.978, 0.923, 0.709, 0.99, 0.958, 0.335, 0.774, 0.84, 0.926, 0.51, 0.704, 0.975, 0.684, 0.226, 0.811, 0.613, 0.715, 0.936, 0.69, 0.987, 0.62, 0.735, 0.989, 0.769, 0.945, 0.926, 0.748, 0.577, 0.585, 0.975, 0.966, 0.357, 0.303, 0.741, 0.995, 0.711, 0.712, 0.929, 0.757, 0.939, 0.91, 0.707, 0.831, 0.637, 0.78, 0.26, 0.441, 0.927, 0.826, 0.088, 0.316, 0.898, 0.572, 0.529, 0.873, 0.515, 0.935, 0.968, 0.905, 0.907, 0.998, 0.522, 0.635, 0.379, 0.977, 0.244, 0.82, 0.986, 0.996, 0.849, 0.721, 0.807, 0.89, 0.526, 0.493, 0.737, 0.995, 0.999, 0.553, 0.998, 0.995, 0.919, 0.961, 0.524, 0.519, 0.304, 0.875, 0.616, 0.773, 0.923, 0.707, 0.971, 0.989, 0.967, 0.623, 0.77, 0.974, 0.691, 0.128, 0.991, 0.3, 0.349, 0.598, 0.902, 0.407, 0.973, 0.879, 0.981, 0.898, 0.965, 0.703, 0.913, 0.357, 0.922, 0.451, 0.997, 0.239, 0.472, 0.969, 0.993, 0.06, 0.668, 0.901, 0.62, 0.745, 0.356, 0.558, 0.989, 0.893, 0.984, 0.836, 0.071, 0.223, 0.969, 0.835, 0.834, 0.974, 0.683, 0.196, 0.808, 0.993, 0.844, 0.5, 0.992, 0.472, 0.949, 0.963, 0.662, 0.449, 0.782, 0.329, 0.665, 0.927, 0.992, 0.363, 0.966, 0.625, 0.96, 0.961, 0.938, 0.565, 0.621, 0.825, 0.836, 0.995, 0.487, 0.52, 0.953, 0.961, 0.992, 0.999, 0.982, 0.931, 0.835, 0.989, 0.982, 0.919, 0.878, 0.697, 0.825, 0.803, 0.601, 0.969, 0.859, 0.382, 0.685, 0.847, 0.765, 0.979, 0.702, 0.829, 0.608, 0.263, 0.881, 0.819, 0.847, 0.868, 0.217, 0.844, 0.893, 0.827, 0.99, 0.97, 0.874, 0.883, 0.982, 0.973, 0.873, 0.745, 0.393, 0.829, 0.922, 0.804, 0.356, 0.596, 0.986, 0.767, 0.891, 0.792, 0.991, 0.916, 0.704, 0.803, 0.76, 0.913, 0.793, 0.914, 0.975, 0.431, 0.35, 0.999, 0.732, 0.801, 0.655, 0.531, 0.618, 0.5, 0.322, 0.662, 0.494, 0.893, 0.972, 0.337, 0.713, 0.262, 0.898, 0.769, 0.811, 0.893, 0.952, 0.945, 0.987]
global origin = 1
global destination = 50