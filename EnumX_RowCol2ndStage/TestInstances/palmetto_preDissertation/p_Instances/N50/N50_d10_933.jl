global arcs = [1 5; 1 10; 2 17; 2 39; 3 4; 3 6; 3 11; 3 38; 3 39; 3 40; 3 46; 3 50; 4 21; 4 24; 4 27; 4 36; 5 8; 5 30; 6 9; 6 19; 6 29; 6 32; 6 34; 7 17; 7 32; 7 40; 8 24; 8 38; 8 45; 8 47; 9 3; 9 5; 9 17; 9 32; 9 34; 9 37; 9 48; 10 11; 10 12; 10 20; 10 27; 10 28; 10 34; 10 36; 10 44; 11 4; 11 17; 11 24; 11 31; 11 38; 11 49; 11 50; 12 3; 12 16; 12 20; 12 25; 12 49; 13 26; 13 42; 14 31; 14 32; 14 42; 15 22; 16 9; 16 18; 17 8; 17 14; 17 26; 17 40; 18 4; 18 29; 18 34; 18 39; 18 41; 19 17; 20 5; 20 15; 20 33; 20 38; 20 44; 21 30; 21 50; 22 14; 22 15; 22 16; 22 27; 22 31; 22 39; 23 6; 23 25; 24 5; 24 26; 24 28; 24 47; 25 4; 25 8; 25 32; 25 33; 25 34; 25 39; 25 42; 25 49; 26 2; 26 9; 26 17; 26 25; 26 30; 26 39; 27 7; 27 30; 27 44; 28 5; 28 34; 28 41; 28 42; 29 22; 29 33; 29 37; 29 39; 29 41; 30 5; 30 9; 30 12; 30 37; 30 38; 31 13; 31 20; 31 24; 31 33; 31 45; 32 3; 32 24; 32 34; 32 46; 33 26; 33 27; 33 43; 33 50; 34 27; 34 33; 35 8; 35 22; 35 23; 35 27; 35 31; 35 44; 36 6; 36 13; 36 16; 36 22; 36 25; 36 34; 36 47; 37 3; 37 4; 37 14; 37 21; 37 24; 37 27; 37 32; 38 19; 38 44; 38 50; 39 5; 39 10; 39 12; 39 16; 39 19; 39 29; 39 50; 40 9; 40 10; 40 13; 40 31; 40 35; 41 7; 41 23; 41 42; 42 8; 42 23; 42 28; 42 29; 42 37; 42 41; 42 46; 42 47; 43 19; 43 21; 43 31; 43 39; 43 47; 44 40; 45 4; 45 5; 45 6; 45 7; 45 14; 45 21; 45 37; 46 27; 46 32; 46 35; 47 4; 47 12; 47 26; 47 33; 47 34; 47 49; 48 4; 48 23; 48 28; 48 32; 48 47; 48 50; 49 2; 49 3; 49 17; 49 29; 49 32; 49 33; 49 40]
global d_x = [8.0, 1.0, 7.0, 8.0, 2.0, 10.0, 8.0, 9.0, 10.0, 6.0, 5.0, 10.0, 6.0, 1.0, 10.0, 7.0, 4.0, 10.0, 4.0, 2.0, 6.0, 2.0, 9.0, 3.0, 1.0, 7.0, 2.0, 10.0, 4.0, 7.0, 8.0, 10.0, 9.0, 10.0, 6.0, 10.0, 3.0, 7.0, 10.0, 8.0, 5.0, 5.0, 6.0, 9.0, 4.0, 4.0, 9.0, 3.0, 10.0, 2.0, 1.0, 10.0, 10.0, 8.0, 2.0, 3.0, 5.0, 7.0, 2.0, 7.0, 6.0, 10.0, 9.0, 4.0, 4.0, 3.0, 3.0, 9.0, 6.0, 1.0, 5.0, 8.0, 8.0, 5.0, 7.0, 6.0, 1.0, 7.0, 9.0, 3.0, 2.0, 8.0, 2.0, 9.0, 9.0, 8.0, 6.0, 4.0, 3.0, 8.0, 8.0, 7.0, 9.0, 9.0, 9.0, 4.0, 10.0, 7.0, 6.0, 2.0, 1.0, 10.0, 5.0, 1.0, 7.0, 9.0, 2.0, 2.0, 3.0, 5.0, 5.0, 4.0, 7.0, 4.0, 4.0, 2.0, 9.0, 7.0, 10.0, 1.0, 10.0, 9.0, 10.0, 3.0, 8.0, 6.0, 4.0, 6.0, 2.0, 6.0, 8.0, 4.0, 3.0, 10.0, 3.0, 5.0, 8.0, 4.0, 2.0, 1.0, 9.0, 8.0, 1.0, 1.0, 4.0, 2.0, 2.0, 1.0, 1.0, 10.0, 4.0, 7.0, 6.0, 6.0, 2.0, 1.0, 9.0, 4.0, 7.0, 4.0, 1.0, 4.0, 2.0, 8.0, 4.0, 5.0, 7.0, 7.0, 4.0, 2.0, 9.0, 6.0, 7.0, 4.0, 3.0, 9.0, 7.0, 2.0, 7.0, 3.0, 5.0, 10.0, 7.0, 3.0, 9.0, 6.0, 7.0, 5.0, 4.0, 8.0, 3.0, 4.0, 1.0, 5.0, 10.0, 1.0, 9.0, 8.0, 2.0, 5.0, 3.0, 10.0, 5.0, 10.0, 2.0, 5.0, 2.0, 6.0, 4.0, 7.0, 9.0, 6.0, 6.0, 10.0, 8.0, 5.0, 9.0, 10.0, 8.0, 8.0, 6.0]
global b_x = 5
global d_y = [5.0, 4.0, 10.0, 8.0, 10.0, 8.0, 3.0, 5.0, 4.0, 2.0, 8.0, 8.0, 8.0, 8.0, 9.0, 1.0, 2.0, 3.0, 3.0, 8.0, 3.0, 5.0, 6.0, 7.0, 3.0, 5.0, 10.0, 6.0, 2.0, 9.0, 3.0, 9.0, 9.0, 5.0, 8.0, 6.0, 8.0, 7.0, 1.0, 1.0, 5.0, 4.0, 8.0, 6.0, 5.0, 9.0, 3.0, 9.0, 5.0, 5.0, 3.0, 7.0, 8.0, 1.0, 7.0, 3.0, 10.0, 4.0, 5.0, 2.0, 7.0, 7.0, 9.0, 1.0, 2.0, 1.0, 6.0, 10.0, 9.0, 8.0, 6.0, 3.0, 6.0, 3.0, 5.0, 9.0, 8.0, 4.0, 8.0, 9.0, 7.0, 5.0, 10.0, 7.0, 8.0, 5.0, 6.0, 8.0, 1.0, 7.0, 7.0, 10.0, 1.0, 3.0, 5.0, 10.0, 8.0, 9.0, 8.0, 10.0, 9.0, 10.0, 8.0, 7.0, 8.0, 4.0, 4.0, 2.0, 7.0, 5.0, 6.0, 6.0, 5.0, 6.0, 8.0, 7.0, 2.0, 1.0, 6.0, 4.0, 3.0, 10.0, 1.0, 1.0, 1.0, 3.0, 10.0, 4.0, 9.0, 9.0, 2.0, 2.0, 5.0, 8.0, 7.0, 10.0, 4.0, 3.0, 1.0, 5.0, 8.0, 3.0, 10.0, 8.0, 4.0, 5.0, 10.0, 2.0, 10.0, 9.0, 3.0, 4.0, 7.0, 4.0, 4.0, 4.0, 8.0, 7.0, 1.0, 10.0, 4.0, 4.0, 7.0, 1.0, 10.0, 9.0, 6.0, 9.0, 2.0, 3.0, 6.0, 2.0, 6.0, 5.0, 5.0, 7.0, 3.0, 6.0, 7.0, 5.0, 3.0, 8.0, 6.0, 6.0, 2.0, 1.0, 4.0, 2.0, 1.0, 6.0, 7.0, 1.0, 4.0, 5.0, 10.0, 5.0, 2.0, 1.0, 7.0, 2.0, 6.0, 4.0, 10.0, 9.0, 2.0, 7.0, 9.0, 3.0, 6.0, 10.0, 9.0, 1.0, 6.0, 8.0, 8.0, 7.0, 9.0, 4.0, 3.0, 2.0, 5.0]
global b_y = 10
global p = [0.001, 0.57, 0.915, 0.383, 0.63, 0.808, 0.813, 0.433, 0.528, 0.801, 0.236, 0.628, 0.38, 0.858, 0.983, 0.252, 0.675, 0.07, 0.164, 0.244, 0.194, 0.572, 0.207, 0.246, 0.556, 0.912, 0.106, 0.133, 0.322, 0.421, 0.731, 0.385, 0.709, 0.17, 0.895, 0.728, 0.664, 0.168, 0.601, 0.229, 0.565, 0.56, 0.377, 0.531, 0.691, 0.791, 0.851, 0.39, 0.438, 0.118, 0.605, 0.725, 0.609, 0.567, 0.888, 0.121, 0.068, 0.44, 0.534, 0.105, 0.193, 0.647, 0.486, 0.776, 0.217, 0.668, 0.601, 0.212, 0.144, 0.401, 0.923, 0.274, 0.713, 0.504, 0.252, 0.251, 0.34, 0.212, 0.998, 0.724, 0.898, 0.82, 0.734, 0.035, 0.178, 0.589, 0.252, 0.332, 0.97, 0.702, 0.223, 0.391, 0.673, 0.418, 0.09, 0.062, 0.999, 0.034, 0.676, 0.182, 0.714, 0.065, 0.538, 0.912, 0.795, 0.854, 0.284, 0.445, 0.823, 0.641, 0.909, 0.439, 0.995, 0.657, 0.635, 0.447, 0.762, 0.842, 0.364, 0.726, 0.764, 0.203, 0.799, 0.949, 0.727, 0.101, 0.357, 0.744, 0.279, 0.788, 0.463, 0.923, 0.844, 0.258, 0.194, 0.389, 0.234, 0.13, 0.388, 0.355, 0.18, 0.221, 0.613, 0.773, 0.202, 0.977, 0.649, 0.062, 0.664, 0.878, 0.626, 0.979, 0.546, 0.987, 0.627, 0.43, 0.386, 0.106, 0.972, 0.05, 0.718, 0.184, 0.948, 0.172, 0.421, 0.491, 0.448, 0.488, 0.821, 0.023, 0.369, 0.315, 0.805, 0.05, 0.122, 0.079, 0.212, 0.174, 0.033, 0.725, 0.037, 0.822, 0.251, 0.692, 0.301, 0.079, 0.646, 0.706, 0.121, 0.936, 0.28, 0.798, 0.85, 0.267, 0.387, 0.687, 0.241, 0.868, 0.344, 0.967, 0.249, 0.334, 0.54, 0.114, 0.017, 0.57, 0.093, 0.901, 0.127, 0.093, 0.745, 0.304, 0.091, 0.077, 0.475, 0.516, 0.425, 0.68, 0.879, 0.362, 0.571]
global q = [0.733, 0.927, 0.989, 0.556, 0.738, 0.983, 0.944, 0.484, 0.944, 0.815, 0.432, 0.954, 0.996, 0.999, 0.986, 0.42, 0.885, 0.373, 0.298, 0.713, 0.273, 0.636, 0.54, 0.339, 0.899, 0.947, 0.929, 0.992, 0.936, 0.703, 0.79, 0.822, 0.946, 0.813, 0.931, 0.991, 0.818, 0.205, 0.759, 0.746, 0.628, 0.998, 0.658, 0.602, 0.728, 0.962, 0.922, 0.668, 0.514, 0.731, 0.943, 0.912, 0.641, 0.709, 0.964, 0.202, 0.865, 0.915, 0.987, 0.839, 0.789, 0.738, 0.714, 0.858, 0.736, 0.95, 0.676, 0.31, 0.804, 0.717, 0.928, 0.711, 0.947, 0.989, 0.709, 0.863, 0.571, 0.431, 0.999, 0.988, 0.927, 0.941, 0.867, 0.672, 0.857, 0.872, 0.535, 0.982, 0.972, 0.778, 0.529, 0.828, 0.747, 0.879, 0.621, 0.108, 0.999, 0.118, 0.679, 0.25, 0.884, 0.612, 0.832, 0.932, 0.906, 0.987, 0.771, 0.513, 0.994, 0.79, 0.911, 0.762, 0.997, 0.989, 0.792, 0.963, 0.903, 0.962, 0.396, 0.962, 0.905, 0.318, 0.984, 0.98, 0.749, 0.878, 0.569, 0.99, 0.766, 0.835, 0.546, 0.971, 0.886, 0.721, 0.755, 0.772, 0.864, 0.581, 0.658, 0.833, 0.653, 0.394, 0.994, 0.997, 0.93, 0.993, 0.693, 0.325, 0.877, 0.953, 0.656, 0.989, 0.754, 0.998, 0.889, 0.717, 0.655, 0.471, 0.977, 0.754, 0.88, 0.408, 0.993, 0.273, 0.74, 0.51, 0.903, 0.999, 0.985, 0.06, 0.832, 0.867, 0.946, 0.41, 0.93, 0.295, 0.558, 0.682, 0.377, 0.87, 0.626, 0.932, 0.927, 0.923, 0.419, 0.091, 0.798, 0.973, 0.745, 0.992, 0.662, 0.986, 0.999, 0.883, 0.519, 0.859, 0.532, 0.878, 0.355, 0.977, 0.849, 0.595, 0.625, 0.95, 0.316, 0.679, 0.408, 0.935, 0.481, 0.468, 0.955, 0.624, 0.155, 0.941, 0.761, 0.647, 0.858, 0.92, 0.964, 0.689, 0.787]
global origin = 1
global destination = 50