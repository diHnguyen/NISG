global arcs = [1 7; 1 13; 1 23; 1 27; 2 8; 2 17; 2 23; 2 27; 2 31; 2 35; 3 11; 3 29; 3 40; 4 6; 4 14; 4 25; 4 26; 4 30; 4 31; 4 32; 5 18; 5 25; 5 40; 6 21; 6 22; 6 29; 6 40; 7 3; 7 11; 7 12; 7 22; 7 37; 8 4; 8 24; 8 25; 9 2; 9 24; 9 25; 9 29; 10 2; 10 3; 10 14; 10 15; 10 25; 10 33; 10 40; 11 4; 11 7; 11 9; 11 14; 11 15; 11 30; 12 25; 12 35; 12 38; 13 4; 13 5; 13 15; 13 24; 13 30; 13 38; 14 38; 15 13; 15 38; 16 7; 16 14; 16 24; 16 35; 17 2; 17 29; 17 35; 18 2; 18 11; 18 40; 19 15; 19 21; 20 15; 20 17; 20 22; 20 23; 20 29; 20 30; 20 35; 20 39; 21 8; 21 13; 21 32; 21 39; 22 14; 22 17; 22 38; 23 5; 23 15; 23 34; 24 5; 24 7; 24 9; 24 18; 24 20; 24 40; 25 6; 25 35; 26 6; 26 20; 26 25; 26 32; 26 36; 27 4; 27 13; 28 2; 28 14; 28 25; 28 31; 28 37; 29 4; 29 6; 29 11; 30 15; 30 25; 30 33; 30 35; 31 2; 31 38; 32 11; 33 2; 33 18; 33 22; 33 36; 34 6; 34 22; 34 31; 34 33; 35 10; 35 19; 35 39; 36 4; 36 5; 36 10; 36 11; 36 15; 36 16; 36 18; 36 25; 36 28; 36 30; 36 31; 36 34; 37 9; 37 19; 37 21; 37 25; 37 29; 37 31; 37 36; 38 10; 38 12; 38 29; 39 37]
global d_x = [1.0, 2.0, 9.0, 6.0, 9.0, 4.0, 2.0, 3.0, 10.0, 2.0, 3.0, 6.0, 6.0, 8.0, 6.0, 9.0, 9.0, 3.0, 9.0, 10.0, 4.0, 4.0, 10.0, 2.0, 8.0, 5.0, 1.0, 9.0, 1.0, 10.0, 4.0, 5.0, 9.0, 1.0, 8.0, 2.0, 6.0, 3.0, 1.0, 6.0, 10.0, 1.0, 7.0, 8.0, 2.0, 3.0, 7.0, 2.0, 1.0, 8.0, 5.0, 4.0, 2.0, 8.0, 4.0, 10.0, 4.0, 2.0, 5.0, 4.0, 1.0, 2.0, 3.0, 4.0, 4.0, 9.0, 1.0, 10.0, 3.0, 7.0, 4.0, 10.0, 9.0, 4.0, 6.0, 5.0, 1.0, 7.0, 8.0, 2.0, 2.0, 7.0, 1.0, 3.0, 5.0, 8.0, 10.0, 2.0, 5.0, 9.0, 3.0, 9.0, 10.0, 8.0, 1.0, 8.0, 3.0, 1.0, 10.0, 2.0, 9.0, 5.0, 4.0, 3.0, 10.0, 5.0, 1.0, 7.0, 4.0, 3.0, 7.0, 4.0, 9.0, 9.0, 4.0, 8.0, 7.0, 9.0, 1.0, 8.0, 2.0, 3.0, 8.0, 8.0, 3.0, 6.0, 1.0, 9.0, 5.0, 6.0, 3.0, 4.0, 10.0, 1.0, 8.0, 4.0, 8.0, 7.0, 3.0, 10.0, 9.0, 7.0, 6.0, 8.0, 3.0, 6.0, 3.0, 6.0, 7.0, 10.0, 8.0, 9.0, 1.0, 3.0, 2.0, 2.0, 6.0, 3.0]
global b_x = 5
global d_y = [6.0, 3.0, 10.0, 4.0, 3.0, 5.0, 1.0, 2.0, 4.0, 2.0, 4.0, 1.0, 7.0, 8.0, 8.0, 3.0, 8.0, 2.0, 8.0, 10.0, 8.0, 4.0, 9.0, 8.0, 8.0, 1.0, 2.0, 3.0, 2.0, 3.0, 10.0, 7.0, 9.0, 4.0, 10.0, 1.0, 8.0, 8.0, 2.0, 3.0, 2.0, 8.0, 3.0, 5.0, 2.0, 1.0, 2.0, 1.0, 6.0, 5.0, 6.0, 6.0, 5.0, 9.0, 3.0, 10.0, 6.0, 3.0, 3.0, 9.0, 9.0, 2.0, 5.0, 8.0, 7.0, 7.0, 6.0, 9.0, 1.0, 5.0, 7.0, 7.0, 5.0, 5.0, 7.0, 8.0, 3.0, 9.0, 8.0, 10.0, 9.0, 7.0, 4.0, 10.0, 10.0, 5.0, 3.0, 8.0, 5.0, 10.0, 7.0, 8.0, 7.0, 7.0, 5.0, 8.0, 2.0, 6.0, 6.0, 3.0, 9.0, 9.0, 1.0, 4.0, 9.0, 10.0, 10.0, 4.0, 2.0, 4.0, 10.0, 5.0, 6.0, 1.0, 8.0, 3.0, 5.0, 9.0, 9.0, 4.0, 5.0, 8.0, 9.0, 8.0, 1.0, 5.0, 2.0, 4.0, 5.0, 6.0, 7.0, 10.0, 3.0, 1.0, 10.0, 4.0, 10.0, 3.0, 3.0, 4.0, 1.0, 3.0, 10.0, 8.0, 9.0, 3.0, 4.0, 2.0, 3.0, 8.0, 9.0, 3.0, 3.0, 9.0, 7.0, 10.0, 2.0, 1.0]
global b_y = 10
global p = [0.593, 0.997, 0.404, 0.914, 0.828, 0.298, 0.47, 0.017, 0.471, 0.348, 0.867, 0.684, 0.24, 0.574, 0.732, 0.002, 0.368, 0.193, 0.506, 0.789, 0.748, 0.307, 0.985, 0.994, 0.269, 0.446, 0.76, 0.59, 0.802, 0.126, 0.921, 0.785, 0.356, 0.828, 0.95, 0.179, 0.206, 0.234, 0.448, 0.27, 0.115, 0.691, 0.422, 0.054, 0.188, 0.657, 0.488, 0.669, 0.391, 0.985, 0.107, 0.073, 0.469, 0.815, 0.293, 0.823, 0.888, 0.774, 0.233, 0.557, 0.204, 0.064, 0.006, 0.417, 0.105, 0.883, 0.397, 0.521, 0.189, 0.694, 0.171, 0.688, 0.443, 0.471, 0.609, 0.642, 0.4, 0.44, 0.007, 0.95, 0.726, 0.229, 0.502, 0.631, 0.548, 0.414, 0.71, 0.433, 0.986, 0.609, 0.856, 0.637, 0.159, 0.445, 0.059, 0.627, 0.622, 0.672, 0.503, 0.439, 0.346, 0.645, 0.703, 0.504, 0.439, 0.683, 0.934, 0.451, 0.349, 0.452, 0.83, 0.437, 0.963, 0.718, 0.223, 0.091, 0.851, 0.169, 0.889, 0.129, 0.299, 0.37, 0.297, 0.185, 0.526, 0.776, 0.742, 0.602, 0.281, 0.395, 0.702, 0.608, 0.351, 0.973, 0.741, 0.029, 0.146, 0.734, 0.383, 0.32, 0.126, 0.609, 0.984, 0.193, 0.298, 0.036, 0.21, 0.891, 0.999, 0.727, 0.571, 0.993, 0.162, 0.931, 0.866, 0.135, 0.622, 0.118]
global q = [0.907, 0.999, 0.703, 0.939, 0.991, 0.664, 0.606, 0.712, 0.738, 0.83, 0.902, 0.943, 0.665, 0.73, 0.844, 0.254, 0.58, 0.833, 0.587, 0.872, 0.984, 0.619, 0.985, 0.996, 0.368, 0.83, 0.861, 0.674, 0.957, 0.843, 0.994, 0.876, 0.361, 0.969, 0.963, 0.388, 0.95, 0.751, 0.607, 0.28, 0.191, 0.77, 0.55, 0.306, 0.254, 0.733, 0.648, 0.721, 0.824, 0.987, 0.164, 0.229, 0.977, 0.871, 0.639, 0.896, 0.902, 0.918, 0.339, 0.772, 0.678, 0.219, 0.828, 0.54, 0.796, 0.898, 0.971, 0.772, 0.928, 0.81, 0.688, 0.984, 0.576, 0.587, 0.692, 0.76, 0.558, 0.919, 0.632, 0.987, 0.846, 0.294, 0.706, 0.964, 0.623, 0.894, 0.79, 0.621, 0.992, 0.613, 0.999, 0.79, 0.215, 0.94, 0.462, 0.654, 0.707, 0.86, 0.645, 0.934, 0.367, 0.933, 0.995, 0.715, 0.562, 0.999, 0.967, 0.65, 0.86, 0.53, 0.921, 0.98, 0.999, 0.928, 0.872, 0.905, 0.952, 0.468, 0.957, 0.56, 0.816, 0.584, 0.376, 0.755, 0.783, 0.935, 0.932, 0.806, 0.586, 0.853, 0.727, 0.936, 0.542, 0.999, 0.764, 0.77, 0.35, 0.812, 0.785, 0.811, 0.689, 0.892, 0.988, 0.596, 0.756, 0.045, 0.256, 0.914, 0.999, 0.775, 0.775, 0.993, 0.907, 0.931, 0.902, 0.784, 0.995, 0.691]
global origin = 1
global destination = 40