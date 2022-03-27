global arcs = [1 4; 1 6; 1 12; 1 16; 1 19; 1 22; 1 26; 1 28; 2 4; 2 6; 2 9; 2 11; 2 16; 2 17; 2 18; 2 24; 3 7; 3 10; 3 12; 3 22; 3 27; 4 3; 4 9; 4 15; 4 16; 4 27; 5 2; 5 4; 5 7; 5 10; 5 19; 5 21; 6 5; 6 9; 6 19; 6 24; 6 29; 6 30; 7 4; 7 5; 7 10; 7 14; 7 17; 8 2; 8 14; 8 18; 9 3; 9 5; 9 6; 9 12; 9 23; 9 26; 9 27; 9 29; 10 6; 10 7; 10 9; 10 13; 10 14; 10 28; 11 5; 11 29; 12 8; 12 11; 12 14; 12 16; 12 18; 12 21; 12 23; 13 8; 13 11; 13 23; 13 24; 13 30; 14 3; 14 8; 14 9; 14 10; 14 15; 14 19; 15 3; 15 11; 15 19; 15 23; 15 26; 16 14; 16 15; 16 21; 16 30; 17 5; 17 9; 17 10; 17 13; 17 14; 17 22; 17 23; 18 6; 18 9; 18 12; 18 20; 18 23; 18 25; 18 26; 18 27; 19 10; 20 4; 20 5; 20 14; 20 18; 20 22; 20 24; 20 25; 21 2; 21 3; 21 4; 21 7; 21 8; 21 11; 21 13; 21 20; 21 27; 21 29; 22 2; 22 5; 22 18; 22 19; 22 25; 23 3; 23 4; 23 8; 23 30; 24 13; 24 20; 24 22; 24 29; 25 14; 25 20; 25 21; 25 23; 26 6; 26 15; 26 17; 26 19; 26 25; 27 5; 27 8; 27 11; 27 12; 27 15; 28 2; 28 5; 28 18; 28 21; 28 27; 29 4; 29 6; 29 13; 29 21; 29 23; 29 26]
global d_x = [10.0, 6.0, 10.0, 3.0, 4.0, 9.0, 6.0, 5.0, 10.0, 6.0, 1.0, 9.0, 6.0, 3.0, 1.0, 8.0, 8.0, 7.0, 6.0, 10.0, 5.0, 1.0, 1.0, 6.0, 5.0, 5.0, 3.0, 1.0, 3.0, 2.0, 1.0, 7.0, 5.0, 2.0, 6.0, 1.0, 2.0, 9.0, 2.0, 6.0, 2.0, 7.0, 2.0, 9.0, 8.0, 7.0, 3.0, 8.0, 10.0, 3.0, 9.0, 2.0, 10.0, 1.0, 1.0, 6.0, 3.0, 1.0, 7.0, 8.0, 6.0, 8.0, 3.0, 5.0, 5.0, 7.0, 3.0, 8.0, 6.0, 7.0, 1.0, 9.0, 1.0, 7.0, 4.0, 3.0, 4.0, 10.0, 2.0, 5.0, 2.0, 3.0, 2.0, 3.0, 3.0, 9.0, 7.0, 1.0, 7.0, 4.0, 6.0, 9.0, 7.0, 10.0, 6.0, 1.0, 4.0, 8.0, 10.0, 7.0, 3.0, 7.0, 10.0, 2.0, 8.0, 4.0, 8.0, 4.0, 10.0, 3.0, 10.0, 5.0, 5.0, 4.0, 5.0, 1.0, 7.0, 6.0, 2.0, 7.0, 3.0, 8.0, 2.0, 6.0, 6.0, 6.0, 10.0, 4.0, 4.0, 8.0, 10.0, 2.0, 2.0, 5.0, 1.0, 8.0, 6.0, 5.0, 2.0, 8.0, 8.0, 9.0, 1.0, 6.0, 5.0, 10.0, 1.0, 2.0, 2.0, 10.0, 10.0, 6.0, 8.0, 7.0, 4.0, 8.0, 4.0, 2.0, 5.0, 8.0]
global b_x = 5
global d_y = [5.0, 7.0, 10.0, 8.0, 8.0, 4.0, 5.0, 6.0, 5.0, 2.0, 6.0, 1.0, 8.0, 2.0, 2.0, 3.0, 6.0, 6.0, 9.0, 3.0, 5.0, 3.0, 5.0, 6.0, 5.0, 7.0, 2.0, 9.0, 3.0, 6.0, 9.0, 8.0, 5.0, 10.0, 8.0, 3.0, 8.0, 9.0, 10.0, 1.0, 1.0, 7.0, 3.0, 10.0, 8.0, 1.0, 10.0, 8.0, 2.0, 9.0, 9.0, 1.0, 4.0, 6.0, 8.0, 4.0, 8.0, 6.0, 9.0, 5.0, 7.0, 7.0, 6.0, 1.0, 1.0, 2.0, 3.0, 10.0, 2.0, 9.0, 4.0, 8.0, 6.0, 3.0, 4.0, 4.0, 9.0, 5.0, 5.0, 9.0, 7.0, 6.0, 2.0, 9.0, 10.0, 1.0, 7.0, 9.0, 3.0, 7.0, 5.0, 1.0, 1.0, 7.0, 1.0, 8.0, 2.0, 9.0, 8.0, 7.0, 9.0, 9.0, 6.0, 6.0, 7.0, 3.0, 7.0, 10.0, 2.0, 10.0, 2.0, 4.0, 8.0, 7.0, 4.0, 7.0, 3.0, 8.0, 7.0, 9.0, 8.0, 4.0, 2.0, 5.0, 9.0, 3.0, 9.0, 2.0, 6.0, 2.0, 4.0, 6.0, 6.0, 5.0, 1.0, 9.0, 4.0, 6.0, 5.0, 10.0, 8.0, 10.0, 5.0, 7.0, 10.0, 3.0, 2.0, 8.0, 10.0, 4.0, 4.0, 8.0, 1.0, 4.0, 5.0, 6.0, 6.0, 7.0, 7.0, 8.0]
global b_y = 10
global p = [0.247, 0.599, 0.641, 0.327, 0.831, 0.312, 0.712, 0.761, 0.919, 0.117, 0.512, 0.015, 0.124, 0.698, 0.124, 0.739, 0.852, 0.222, 0.969, 0.57, 0.085, 0.458, 0.96, 0.93, 0.415, 0.442, 0.407, 0.823, 0.95, 0.479, 0.043, 0.877, 0.095, 0.767, 0.446, 0.569, 0.364, 0.743, 0.413, 0.269, 0.899, 0.749, 0.649, 0.978, 0.393, 0.164, 0.509, 0.786, 0.055, 0.59, 0.491, 0.876, 0.197, 0.359, 0.762, 0.432, 0.795, 0.456, 0.646, 0.812, 0.312, 0.237, 0.618, 0.24, 0.156, 0.538, 0.027, 0.344, 0.382, 0.507, 0.237, 0.696, 0.842, 0.523, 0.815, 0.671, 0.083, 0.145, 0.009, 0.696, 0.86, 0.625, 0.284, 0.951, 0.752, 0.407, 0.795, 0.629, 0.24, 0.967, 0.993, 0.545, 0.461, 0.602, 0.161, 0.901, 0.105, 0.791, 0.145, 0.585, 0.435, 0.335, 0.015, 0.591, 0.113, 0.322, 0.66, 0.772, 0.963, 0.699, 0.835, 0.36, 0.491, 0.99, 0.893, 0.751, 0.49, 0.155, 0.811, 0.366, 0.931, 0.292, 0.562, 0.222, 0.202, 0.212, 0.98, 0.587, 0.628, 0.605, 0.914, 0.885, 0.705, 0.634, 0.083, 0.994, 0.44, 0.128, 0.715, 0.679, 0.886, 0.598, 0.758, 0.439, 0.023, 0.75, 0.636, 0.543, 0.044, 0.79, 0.23, 0.032, 0.021, 0.073, 0.846, 0.277, 0.469, 0.271, 0.54, 0.362]
global q = [0.622, 0.701, 0.819, 0.465, 0.925, 0.318, 0.762, 0.893, 0.987, 0.196, 0.741, 0.461, 0.526, 0.88, 0.621, 0.853, 0.997, 0.732, 0.984, 0.942, 0.93, 0.633, 0.987, 0.977, 0.96, 0.705, 0.661, 0.974, 0.994, 0.948, 0.562, 0.93, 0.887, 0.962, 0.942, 0.733, 0.579, 0.754, 0.467, 0.42, 0.978, 0.947, 0.979, 0.981, 0.475, 0.614, 0.732, 0.885, 0.546, 0.96, 0.909, 0.92, 0.661, 0.658, 0.968, 0.946, 0.911, 0.715, 0.855, 0.894, 0.435, 0.318, 0.695, 0.825, 0.513, 0.969, 0.214, 0.385, 0.981, 0.58, 0.481, 0.9, 0.937, 0.921, 0.819, 0.904, 0.859, 0.623, 0.126, 0.86, 0.901, 0.992, 0.366, 0.975, 0.786, 0.758, 0.933, 0.7, 0.273, 0.997, 0.999, 0.724, 0.953, 0.706, 0.641, 0.917, 0.874, 0.911, 0.9, 0.939, 0.964, 0.525, 0.336, 0.855, 0.31, 0.974, 0.785, 0.834, 0.977, 0.82, 0.855, 0.906, 0.891, 0.998, 0.983, 0.817, 0.578, 0.451, 0.838, 0.473, 0.977, 0.499, 0.772, 0.76, 0.208, 0.379, 0.981, 0.942, 0.713, 0.865, 0.947, 0.952, 0.975, 0.906, 0.882, 0.997, 0.867, 0.501, 0.763, 0.681, 0.992, 0.998, 0.789, 0.485, 0.516, 0.999, 0.731, 0.596, 0.949, 0.957, 0.805, 0.619, 0.351, 0.46, 0.853, 0.911, 0.637, 0.88, 0.85, 0.989]
global origin = 1
global destination = 30