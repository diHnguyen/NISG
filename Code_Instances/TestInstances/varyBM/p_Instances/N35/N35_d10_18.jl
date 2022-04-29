global arcs = [1 4; 1 22; 1 25; 1 28; 2 3; 2 4; 2 7; 2 21; 2 28; 2 30; 3 15; 3 29; 3 31; 4 2; 4 7; 4 8; 4 14; 4 19; 4 30; 4 33; 4 35; 5 4; 5 16; 5 19; 5 25; 6 5; 6 16; 6 18; 6 25; 7 10; 7 28; 7 29; 8 3; 8 25; 9 5; 9 20; 9 24; 9 26; 9 35; 10 5; 10 15; 10 31; 11 17; 11 28; 12 26; 13 5; 13 14; 13 22; 13 33; 13 35; 14 32; 15 2; 15 17; 16 8; 16 12; 16 19; 16 20; 17 4; 17 7; 17 13; 17 14; 17 15; 17 16; 17 19; 17 21; 17 23; 18 15; 18 29; 18 34; 19 7; 19 28; 20 3; 20 28; 20 35; 21 7; 21 30; 22 8; 22 12; 23 9; 23 16; 23 31; 24 10; 24 13; 24 16; 24 23; 24 29; 25 3; 25 21; 25 22; 25 27; 25 35; 26 5; 26 6; 26 11; 26 19; 26 20; 26 21; 26 32; 27 14; 28 5; 28 9; 28 13; 28 35; 29 16; 29 30; 29 33; 30 15; 30 25; 30 26; 31 5; 31 17; 31 30; 31 32; 32 5; 32 7; 32 28; 32 30; 33 2; 33 3; 33 22; 33 24; 33 25; 33 27; 33 31; 34 9; 34 15; 34 22; 34 30; 34 32]
global d_x = [4.0, 1.0, 5.0, 7.0, 1.0, 2.0, 6.0, 9.0, 6.0, 3.0, 2.0, 1.0, 4.0, 2.0, 5.0, 4.0, 3.0, 8.0, 3.0, 8.0, 7.0, 1.0, 1.0, 6.0, 2.0, 5.0, 10.0, 6.0, 4.0, 10.0, 1.0, 10.0, 4.0, 3.0, 10.0, 6.0, 1.0, 8.0, 2.0, 5.0, 3.0, 5.0, 6.0, 7.0, 4.0, 6.0, 2.0, 6.0, 10.0, 6.0, 5.0, 2.0, 3.0, 4.0, 1.0, 7.0, 9.0, 10.0, 9.0, 1.0, 7.0, 1.0, 5.0, 2.0, 9.0, 8.0, 3.0, 1.0, 2.0, 1.0, 8.0, 3.0, 1.0, 3.0, 5.0, 8.0, 9.0, 6.0, 2.0, 6.0, 4.0, 5.0, 8.0, 7.0, 5.0, 2.0, 9.0, 3.0, 1.0, 5.0, 10.0, 7.0, 1.0, 6.0, 4.0, 10.0, 4.0, 1.0, 7.0, 8.0, 8.0, 8.0, 5.0, 10.0, 8.0, 7.0, 2.0, 10.0, 10.0, 10.0, 1.0, 2.0, 10.0, 2.0, 2.0, 1.0, 2.0, 10.0, 3.0, 9.0, 10.0, 9.0, 4.0, 6.0, 3.0, 10.0, 9.0, 6.0, 2.0]
global b_x = 5
global d_y = [3.0, 6.0, 10.0, 8.0, 10.0, 5.0, 10.0, 3.0, 4.0, 9.0, 10.0, 2.0, 6.0, 4.0, 10.0, 4.0, 9.0, 2.0, 2.0, 4.0, 2.0, 3.0, 2.0, 5.0, 2.0, 7.0, 1.0, 2.0, 1.0, 1.0, 2.0, 6.0, 8.0, 2.0, 7.0, 7.0, 3.0, 8.0, 3.0, 7.0, 9.0, 8.0, 1.0, 9.0, 1.0, 5.0, 3.0, 2.0, 6.0, 5.0, 1.0, 4.0, 6.0, 1.0, 10.0, 5.0, 9.0, 7.0, 10.0, 5.0, 8.0, 4.0, 6.0, 6.0, 7.0, 2.0, 8.0, 9.0, 3.0, 7.0, 7.0, 7.0, 2.0, 2.0, 3.0, 8.0, 1.0, 2.0, 4.0, 3.0, 6.0, 9.0, 1.0, 4.0, 7.0, 8.0, 2.0, 7.0, 2.0, 9.0, 10.0, 10.0, 2.0, 5.0, 6.0, 5.0, 10.0, 9.0, 4.0, 6.0, 8.0, 8.0, 10.0, 1.0, 2.0, 2.0, 6.0, 10.0, 3.0, 4.0, 9.0, 9.0, 2.0, 3.0, 6.0, 10.0, 5.0, 5.0, 6.0, 5.0, 8.0, 4.0, 6.0, 10.0, 10.0, 9.0, 1.0, 1.0, 6.0]
global b_y = 10
global p = [0.55, 0.715, 0.88, 0.805, 0.985, 0.279, 0.85, 0.231, 0.846, 0.919, 0.727, 0.135, 0.204, 0.356, 0.135, 0.099, 0.816, 0.789, 0.54, 0.013, 0.886, 0.589, 0.146, 0.997, 0.373, 0.033, 0.901, 0.947, 0.156, 0.408, 0.771, 0.766, 0.247, 0.47, 0.746, 0.98, 0.079, 0.958, 0.232, 0.535, 0.312, 0.057, 0.766, 0.796, 0.032, 0.332, 0.693, 0.819, 0.877, 0.69, 0.456, 0.131, 0.417, 0.302, 0.301, 0.907, 0.601, 0.091, 0.873, 0.823, 0.71, 0.123, 0.586, 0.87, 0.972, 0.188, 0.393, 0.638, 0.456, 0.23, 0.095, 0.33, 0.805, 0.866, 0.805, 0.997, 0.118, 0.492, 0.28, 0.813, 0.155, 0.104, 0.811, 0.819, 0.042, 0.85, 0.509, 0.806, 0.346, 0.954, 0.164, 0.315, 0.654, 0.912, 0.938, 0.994, 0.002, 0.029, 0.239, 0.345, 0.082, 0.742, 0.706, 0.517, 0.416, 0.913, 0.351, 0.424, 0.742, 0.191, 0.346, 0.624, 0.347, 0.411, 0.1, 0.88, 0.2, 0.428, 0.666, 0.392, 0.391, 0.468, 0.266, 0.216, 0.018, 0.123, 0.527, 0.111, 0.112]
global q = [0.85, 0.99, 0.953, 0.889, 0.999, 0.426, 0.856, 0.288, 0.959, 0.983, 0.999, 0.669, 0.434, 0.636, 0.748, 0.525, 0.972, 0.917, 0.96, 0.671, 0.958, 0.879, 0.875, 0.997, 0.641, 0.196, 0.998, 0.952, 0.217, 0.684, 0.948, 0.945, 0.768, 0.521, 0.849, 0.993, 0.858, 0.969, 0.531, 0.679, 0.974, 0.68, 0.948, 0.949, 0.37, 0.37, 0.895, 0.862, 0.889, 0.881, 0.933, 0.697, 0.646, 0.669, 0.42, 0.966, 0.712, 0.393, 0.895, 0.905, 0.93, 0.835, 0.776, 0.96, 0.983, 0.223, 0.806, 0.776, 0.984, 0.305, 0.351, 0.864, 0.814, 0.874, 0.957, 0.997, 0.914, 0.906, 0.902, 0.881, 0.218, 0.321, 0.934, 0.893, 0.821, 0.915, 0.897, 0.933, 0.464, 0.981, 0.81, 0.633, 0.884, 0.964, 0.995, 0.998, 0.889, 0.43, 0.699, 0.728, 0.424, 0.918, 0.922, 0.986, 0.626, 0.991, 0.477, 0.506, 0.796, 0.517, 0.802, 0.624, 0.52, 0.697, 0.807, 0.997, 0.509, 0.844, 0.68, 0.567, 0.965, 0.858, 0.7, 0.511, 0.653, 0.277, 0.846, 0.785, 0.819]
global origin = 1
global destination = 35
