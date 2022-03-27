global arcs = [1 2; 1 4; 1 31; 1 48; 2 10; 2 18; 2 38; 2 40; 3 12; 3 15; 3 31; 3 41; 3 43; 4 6; 4 26; 4 32; 4 40; 5 12; 5 19; 5 36; 5 46; 6 37; 7 2; 7 12; 7 19; 7 37; 7 40; 7 49; 8 7; 8 10; 8 16; 8 31; 8 33; 8 43; 8 44; 8 46; 9 7; 9 11; 9 47; 10 8; 10 12; 10 19; 10 33; 10 41; 10 49; 11 9; 11 24; 11 33; 11 37; 11 38; 11 47; 11 48; 12 8; 12 31; 12 32; 12 36; 12 42; 13 3; 13 9; 13 25; 13 36; 13 45; 13 47; 14 5; 14 20; 14 37; 15 28; 15 35; 15 48; 16 2; 16 6; 16 27; 16 41; 17 15; 17 24; 17 46; 17 48; 17 49; 18 2; 18 3; 18 8; 18 11; 18 23; 18 30; 18 33; 19 4; 19 17; 19 31; 19 49; 20 16; 20 18; 20 36; 20 37; 20 40; 21 3; 21 9; 21 12; 21 13; 21 17; 21 23; 21 25; 21 31; 21 38; 21 39; 22 4; 22 11; 22 18; 22 24; 23 13; 23 18; 23 19; 23 31; 23 43; 23 48; 24 2; 24 14; 24 15; 24 32; 25 14; 25 29; 26 15; 26 16; 26 20; 27 16; 27 26; 27 35; 27 50; 28 43; 29 2; 29 3; 29 8; 29 11; 30 18; 30 21; 30 26; 30 31; 30 45; 31 4; 31 15; 31 30; 31 36; 31 42; 32 2; 32 7; 32 8; 33 8; 33 10; 33 11; 33 12; 33 13; 33 31; 34 3; 34 25; 34 35; 34 38; 34 45; 35 4; 35 12; 35 23; 35 31; 35 37; 36 7; 36 30; 36 42; 36 50; 37 39; 37 40; 38 2; 38 23; 38 34; 39 34; 39 41; 40 11; 40 12; 40 25; 41 3; 41 6; 41 7; 41 19; 41 22; 41 34; 41 39; 42 7; 42 8; 42 17; 42 19; 42 20; 42 50; 43 11; 43 37; 43 45; 43 48; 44 8; 44 11; 44 15; 44 25; 44 37; 44 42; 45 5; 45 14; 45 38; 45 49; 46 10; 46 11; 46 12; 46 25; 46 50; 47 2; 47 11; 47 21; 47 34; 47 48; 48 8; 48 17; 48 19; 48 29; 48 44; 49 15; 49 21; 49 22; 49 26; 49 29; 49 50]
global d_x = [8.0, 2.0, 7.0, 1.0, 6.0, 1.0, 2.0, 9.0, 9.0, 8.0, 5.0, 6.0, 2.0, 7.0, 4.0, 9.0, 1.0, 8.0, 10.0, 10.0, 4.0, 8.0, 5.0, 10.0, 9.0, 6.0, 2.0, 5.0, 6.0, 8.0, 1.0, 10.0, 5.0, 6.0, 9.0, 7.0, 2.0, 7.0, 4.0, 1.0, 6.0, 8.0, 6.0, 3.0, 4.0, 7.0, 4.0, 3.0, 8.0, 2.0, 3.0, 2.0, 5.0, 7.0, 1.0, 6.0, 1.0, 6.0, 10.0, 4.0, 4.0, 7.0, 3.0, 6.0, 10.0, 5.0, 6.0, 2.0, 2.0, 9.0, 4.0, 3.0, 2.0, 5.0, 4.0, 5.0, 3.0, 2.0, 4.0, 8.0, 5.0, 2.0, 4.0, 2.0, 1.0, 2.0, 3.0, 2.0, 6.0, 1.0, 8.0, 8.0, 6.0, 10.0, 7.0, 4.0, 7.0, 10.0, 5.0, 9.0, 8.0, 9.0, 10.0, 4.0, 2.0, 1.0, 1.0, 10.0, 6.0, 3.0, 7.0, 7.0, 2.0, 8.0, 6.0, 7.0, 1.0, 4.0, 10.0, 9.0, 8.0, 4.0, 7.0, 5.0, 5.0, 6.0, 9.0, 10.0, 5.0, 3.0, 7.0, 6.0, 6.0, 1.0, 9.0, 9.0, 5.0, 10.0, 2.0, 8.0, 2.0, 7.0, 4.0, 2.0, 10.0, 4.0, 10.0, 9.0, 1.0, 4.0, 7.0, 5.0, 7.0, 10.0, 3.0, 3.0, 9.0, 7.0, 4.0, 10.0, 2.0, 10.0, 10.0, 4.0, 4.0, 1.0, 1.0, 7.0, 2.0, 2.0, 7.0, 8.0, 3.0, 10.0, 4.0, 1.0, 4.0, 5.0, 2.0, 4.0, 9.0, 2.0, 4.0, 9.0, 2.0, 8.0, 4.0, 4.0, 10.0, 1.0, 10.0, 3.0, 3.0, 3.0, 6.0, 6.0, 7.0, 9.0, 5.0, 6.0, 7.0, 1.0, 8.0, 5.0, 10.0, 6.0, 3.0, 8.0, 4.0, 3.0, 8.0, 8.0, 3.0, 4.0, 1.0, 1.0, 7.0, 9.0, 2.0, 1.0, 6.0, 4.0, 4.0]
global b_x = 5
global d_y = [10.0, 4.0, 8.0, 4.0, 3.0, 10.0, 2.0, 9.0, 6.0, 5.0, 3.0, 10.0, 5.0, 8.0, 5.0, 8.0, 7.0, 3.0, 2.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 10.0, 1.0, 4.0, 4.0, 1.0, 9.0, 4.0, 4.0, 2.0, 1.0, 2.0, 8.0, 5.0, 6.0, 9.0, 1.0, 10.0, 1.0, 3.0, 5.0, 6.0, 8.0, 8.0, 7.0, 2.0, 9.0, 3.0, 4.0, 2.0, 1.0, 9.0, 2.0, 1.0, 9.0, 4.0, 4.0, 3.0, 2.0, 7.0, 6.0, 5.0, 3.0, 7.0, 7.0, 9.0, 10.0, 9.0, 8.0, 2.0, 8.0, 7.0, 6.0, 7.0, 8.0, 9.0, 1.0, 5.0, 6.0, 7.0, 6.0, 6.0, 1.0, 1.0, 1.0, 10.0, 6.0, 8.0, 10.0, 1.0, 8.0, 6.0, 10.0, 2.0, 3.0, 6.0, 2.0, 8.0, 6.0, 4.0, 9.0, 4.0, 10.0, 10.0, 3.0, 4.0, 10.0, 7.0, 6.0, 4.0, 7.0, 10.0, 10.0, 9.0, 4.0, 10.0, 4.0, 10.0, 6.0, 7.0, 8.0, 3.0, 5.0, 7.0, 3.0, 7.0, 3.0, 7.0, 1.0, 1.0, 4.0, 6.0, 2.0, 7.0, 5.0, 7.0, 7.0, 10.0, 2.0, 5.0, 2.0, 4.0, 10.0, 3.0, 8.0, 9.0, 4.0, 1.0, 4.0, 7.0, 8.0, 3.0, 8.0, 4.0, 4.0, 6.0, 1.0, 10.0, 4.0, 10.0, 6.0, 8.0, 10.0, 1.0, 9.0, 1.0, 2.0, 8.0, 6.0, 1.0, 8.0, 7.0, 5.0, 8.0, 7.0, 5.0, 8.0, 6.0, 6.0, 8.0, 4.0, 6.0, 5.0, 4.0, 2.0, 5.0, 3.0, 3.0, 9.0, 3.0, 1.0, 3.0, 8.0, 5.0, 3.0, 6.0, 9.0, 10.0, 8.0, 7.0, 6.0, 2.0, 4.0, 7.0, 6.0, 3.0, 1.0, 7.0, 2.0, 6.0, 10.0, 1.0, 9.0, 3.0, 5.0, 7.0, 3.0, 5.0, 3.0]
global b_y = 10
global p = [0.289, 0.85, 0.892, 0.296, 0.542, 0.547, 0.648, 0.946, 0.172, 0.965, 0.08, 0.451, 0.745, 0.742, 0.375, 0.007, 0.758, 0.257, 0.478, 0.188, 0.461, 0.145, 0.213, 0.769, 0.128, 0.861, 0.114, 0.012, 0.862, 0.957, 0.923, 0.57, 0.248, 0.288, 0.191, 0.043, 0.373, 0.085, 0.809, 0.621, 0.722, 0.079, 0.28, 0.894, 0.227, 0.652, 0.169, 0.788, 0.662, 0.22, 0.133, 0.478, 0.345, 0.657, 0.835, 0.241, 0.074, 0.433, 0.261, 0.969, 0.789, 0.534, 0.03, 0.544, 0.419, 0.579, 0.619, 0.041, 0.78, 0.551, 0.453, 0.604, 0.692, 0.49, 0.664, 0.396, 0.761, 0.293, 0.831, 0.03, 0.537, 0.61, 0.729, 0.756, 0.755, 0.645, 0.402, 0.097, 0.537, 0.692, 0.912, 0.088, 0.114, 0.473, 0.242, 0.938, 0.607, 0.245, 0.391, 0.409, 0.544, 0.45, 0.891, 0.534, 0.943, 0.7, 0.817, 0.719, 0.844, 0.891, 0.755, 0.432, 0.677, 0.144, 0.199, 0.09, 0.941, 0.09, 0.448, 0.946, 0.597, 0.611, 0.778, 0.588, 0.156, 0.115, 0.542, 0.108, 0.663, 0.884, 0.85, 0.983, 0.576, 0.81, 0.306, 0.606, 0.83, 0.813, 0.799, 0.542, 0.205, 0.965, 0.551, 0.417, 0.827, 0.901, 0.959, 0.046, 0.92, 0.892, 0.424, 0.386, 0.045, 0.079, 0.888, 0.589, 0.826, 0.071, 0.519, 0.523, 0.2, 0.967, 0.537, 0.71, 0.956, 0.743, 0.285, 0.096, 0.218, 0.67, 0.444, 0.257, 0.579, 0.189, 0.349, 0.317, 0.129, 0.201, 0.679, 0.515, 0.973, 0.714, 0.078, 0.925, 0.227, 0.536, 0.723, 0.2, 0.668, 0.307, 0.854, 0.748, 0.605, 0.316, 0.81, 0.654, 0.041, 0.536, 0.027, 0.389, 0.567, 0.183, 0.393, 0.418, 0.407, 0.352, 0.024, 0.475, 0.149, 0.095, 0.277, 0.216, 0.17, 0.421, 0.462, 0.819, 0.637, 0.955, 0.979, 0.516, 0.257, 0.275, 0.544]
global q = [0.498, 0.926, 0.919, 0.746, 0.663, 0.813, 0.993, 0.984, 0.304, 0.981, 0.258, 0.918, 0.963, 0.924, 0.424, 0.853, 0.954, 0.268, 0.947, 0.318, 0.601, 0.187, 0.732, 0.986, 0.544, 0.976, 0.97, 0.33, 0.873, 0.974, 0.979, 0.694, 0.357, 0.635, 0.281, 0.323, 0.576, 0.249, 0.954, 0.793, 0.895, 0.943, 0.47, 0.986, 0.675, 0.952, 0.917, 0.901, 0.855, 0.738, 0.156, 0.861, 0.781, 0.676, 0.938, 0.909, 0.742, 0.44, 0.35, 0.991, 0.935, 0.743, 0.849, 0.768, 0.624, 0.975, 0.674, 0.82, 0.985, 0.637, 0.962, 0.876, 0.822, 0.886, 0.913, 0.462, 0.835, 0.517, 0.888, 0.979, 0.682, 0.954, 0.748, 0.816, 0.914, 0.878, 0.819, 0.75, 0.701, 0.7, 0.935, 0.748, 0.664, 0.977, 0.687, 0.947, 0.663, 0.53, 0.795, 0.645, 0.794, 0.776, 0.91, 0.956, 0.994, 0.911, 0.824, 0.938, 0.902, 0.956, 0.936, 0.821, 0.96, 0.74, 0.976, 0.745, 0.997, 0.724, 0.6, 0.996, 0.6, 0.81, 0.933, 0.686, 0.509, 0.127, 0.903, 0.462, 0.985, 0.948, 0.884, 0.993, 0.698, 0.813, 0.481, 0.835, 0.93, 0.868, 0.866, 0.713, 0.446, 0.985, 0.783, 0.587, 0.91, 0.972, 0.982, 0.882, 0.979, 0.985, 0.885, 0.701, 0.101, 0.389, 0.927, 0.697, 0.893, 0.223, 0.844, 0.82, 0.632, 0.999, 0.732, 0.929, 0.965, 0.754, 0.622, 0.126, 0.868, 0.687, 0.569, 0.261, 0.938, 0.33, 0.654, 0.503, 0.869, 0.544, 0.813, 0.793, 0.973, 0.855, 0.82, 0.94, 0.763, 0.818, 0.939, 0.476, 0.86, 0.899, 0.917, 0.788, 0.859, 0.793, 0.932, 0.997, 0.391, 0.669, 0.647, 0.709, 0.728, 0.91, 0.493, 0.481, 0.658, 0.389, 0.506, 0.806, 0.693, 0.397, 0.681, 0.543, 0.725, 0.986, 0.73, 0.879, 0.665, 0.99, 0.987, 0.854, 0.794, 0.475, 0.975]
global origin = 1
global destination = 50