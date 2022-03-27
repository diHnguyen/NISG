global arcs = [1 14; 1 16; 1 29; 2 21; 3 5; 3 8; 3 14; 3 24; 3 27; 3 50; 4 11; 4 19; 4 35; 4 37; 4 44; 5 3; 5 8; 5 19; 5 23; 5 32; 5 33; 5 35; 5 46; 5 49; 6 21; 6 23; 6 26; 6 27; 6 35; 6 36; 7 3; 7 6; 7 10; 7 44; 7 47; 8 7; 8 11; 8 15; 9 7; 9 16; 9 29; 9 40; 9 49; 10 3; 10 14; 10 25; 10 26; 10 31; 10 38; 10 41; 10 46; 11 6; 11 9; 11 13; 11 15; 11 30; 11 42; 11 44; 12 2; 12 4; 12 6; 12 9; 12 22; 12 32; 12 33; 12 35; 12 39; 12 43; 13 8; 13 11; 13 25; 13 26; 13 29; 13 37; 13 41; 13 44; 13 45; 14 26; 15 22; 15 38; 16 12; 16 17; 16 22; 16 28; 17 8; 17 15; 17 30; 17 40; 17 42; 18 25; 18 42; 18 43; 18 48; 19 12; 19 15; 19 25; 19 29; 19 36; 19 39; 19 46; 20 16; 20 46; 20 50; 21 7; 21 26; 21 36; 21 39; 22 9; 22 39; 23 5; 23 8; 23 9; 23 16; 23 26; 24 30; 25 2; 25 6; 25 8; 25 11; 26 22; 26 25; 26 41; 26 42; 26 43; 27 10; 27 14; 27 22; 27 26; 27 31; 27 37; 27 41; 28 6; 28 24; 28 35; 28 46; 29 2; 29 3; 29 10; 29 19; 29 23; 29 27; 29 36; 30 23; 30 35; 30 43; 30 45; 31 2; 31 14; 31 20; 31 21; 31 25; 31 26; 31 44; 31 47; 32 4; 32 10; 32 12; 32 47; 33 5; 33 15; 33 20; 33 28; 33 36; 33 47; 34 7; 34 13; 35 10; 35 24; 35 40; 35 47; 36 23; 37 10; 37 32; 37 39; 37 40; 38 2; 38 21; 38 32; 38 47; 39 47; 40 10; 40 27; 40 29; 40 34; 40 35; 40 50; 41 7; 41 10; 41 18; 42 18; 42 35; 42 48; 43 13; 43 45; 43 46; 43 50; 44 22; 44 40; 45 9; 45 11; 45 12; 45 13; 45 15; 45 20; 45 22; 45 29; 45 34; 46 5; 46 6; 46 20; 46 30; 46 50; 47 3; 47 12; 47 26; 47 28; 47 30; 48 19; 48 22; 48 26; 48 30; 49 7; 49 32; 49 48]
global d_x = [7.0, 7.0, 10.0, 8.0, 3.0, 3.0, 1.0, 4.0, 1.0, 6.0, 8.0, 4.0, 7.0, 3.0, 7.0, 5.0, 6.0, 4.0, 10.0, 4.0, 8.0, 7.0, 10.0, 10.0, 4.0, 6.0, 5.0, 3.0, 3.0, 3.0, 1.0, 8.0, 3.0, 5.0, 9.0, 5.0, 8.0, 8.0, 8.0, 3.0, 5.0, 10.0, 5.0, 6.0, 2.0, 6.0, 10.0, 6.0, 2.0, 7.0, 3.0, 2.0, 8.0, 3.0, 3.0, 3.0, 9.0, 9.0, 2.0, 2.0, 8.0, 7.0, 7.0, 9.0, 9.0, 1.0, 10.0, 7.0, 10.0, 3.0, 8.0, 10.0, 8.0, 4.0, 5.0, 10.0, 10.0, 3.0, 6.0, 4.0, 4.0, 3.0, 10.0, 5.0, 3.0, 5.0, 9.0, 10.0, 7.0, 2.0, 4.0, 4.0, 1.0, 1.0, 4.0, 4.0, 1.0, 9.0, 7.0, 4.0, 2.0, 5.0, 4.0, 1.0, 7.0, 5.0, 4.0, 9.0, 9.0, 2.0, 2.0, 7.0, 1.0, 8.0, 5.0, 9.0, 9.0, 8.0, 2.0, 3.0, 1.0, 4.0, 1.0, 6.0, 7.0, 1.0, 4.0, 10.0, 7.0, 2.0, 8.0, 4.0, 10.0, 3.0, 6.0, 8.0, 8.0, 5.0, 9.0, 8.0, 9.0, 8.0, 10.0, 8.0, 1.0, 3.0, 9.0, 7.0, 6.0, 1.0, 8.0, 9.0, 3.0, 6.0, 7.0, 7.0, 9.0, 1.0, 4.0, 1.0, 2.0, 8.0, 6.0, 7.0, 1.0, 3.0, 9.0, 4.0, 9.0, 2.0, 7.0, 3.0, 10.0, 7.0, 5.0, 6.0, 4.0, 4.0, 5.0, 3.0, 10.0, 8.0, 2.0, 2.0, 10.0, 8.0, 4.0, 5.0, 5.0, 3.0, 3.0, 10.0, 4.0, 5.0, 3.0, 9.0, 1.0, 9.0, 6.0, 10.0, 6.0, 10.0, 8.0, 4.0, 7.0, 7.0, 10.0, 8.0, 6.0, 9.0, 8.0, 5.0, 7.0, 9.0, 1.0, 6.0, 1.0, 8.0, 10.0, 5.0, 9.0, 4.0, 2.0, 10.0]
global b_x = 5
global d_y = [6.0, 1.0, 6.0, 7.0, 6.0, 9.0, 3.0, 2.0, 5.0, 2.0, 3.0, 2.0, 10.0, 6.0, 3.0, 5.0, 3.0, 1.0, 2.0, 8.0, 2.0, 7.0, 2.0, 5.0, 3.0, 1.0, 1.0, 5.0, 4.0, 6.0, 3.0, 10.0, 10.0, 2.0, 2.0, 5.0, 3.0, 4.0, 5.0, 2.0, 5.0, 2.0, 10.0, 1.0, 1.0, 10.0, 1.0, 2.0, 10.0, 6.0, 6.0, 3.0, 2.0, 7.0, 3.0, 2.0, 7.0, 2.0, 5.0, 1.0, 8.0, 4.0, 7.0, 7.0, 1.0, 10.0, 6.0, 9.0, 7.0, 4.0, 9.0, 10.0, 8.0, 10.0, 5.0, 9.0, 6.0, 10.0, 1.0, 5.0, 7.0, 1.0, 8.0, 4.0, 10.0, 8.0, 7.0, 10.0, 3.0, 7.0, 5.0, 4.0, 5.0, 9.0, 1.0, 10.0, 2.0, 3.0, 9.0, 1.0, 4.0, 4.0, 2.0, 9.0, 6.0, 9.0, 3.0, 2.0, 2.0, 7.0, 7.0, 9.0, 7.0, 7.0, 10.0, 4.0, 8.0, 1.0, 1.0, 10.0, 2.0, 7.0, 8.0, 5.0, 8.0, 7.0, 3.0, 6.0, 5.0, 9.0, 3.0, 8.0, 7.0, 3.0, 1.0, 1.0, 6.0, 8.0, 7.0, 7.0, 8.0, 3.0, 8.0, 5.0, 6.0, 4.0, 10.0, 8.0, 9.0, 10.0, 7.0, 8.0, 5.0, 2.0, 1.0, 3.0, 7.0, 1.0, 1.0, 9.0, 9.0, 6.0, 7.0, 6.0, 4.0, 7.0, 1.0, 5.0, 1.0, 4.0, 10.0, 5.0, 1.0, 4.0, 2.0, 10.0, 8.0, 2.0, 6.0, 2.0, 8.0, 3.0, 2.0, 6.0, 3.0, 1.0, 4.0, 7.0, 7.0, 8.0, 6.0, 2.0, 3.0, 7.0, 4.0, 6.0, 3.0, 8.0, 5.0, 3.0, 10.0, 3.0, 9.0, 2.0, 3.0, 3.0, 8.0, 3.0, 6.0, 8.0, 10.0, 3.0, 5.0, 3.0, 3.0, 4.0, 3.0, 5.0, 5.0, 2.0, 2.0, 6.0, 2.0, 2.0]
global b_y = 10
global p = [0.254, 0.92, 0.469, 0.729, 0.65, 0.69, 0.842, 0.561, 0.187, 0.722, 0.632, 0.47, 0.604, 0.729, 0.653, 0.698, 0.072, 0.306, 0.527, 0.346, 0.029, 0.678, 0.19, 0.994, 0.119, 0.232, 0.912, 0.972, 0.952, 0.363, 0.119, 0.564, 0.61, 0.026, 0.862, 0.256, 0.481, 0.619, 0.194, 0.879, 0.446, 0.369, 0.543, 0.21, 0.192, 0.468, 0.978, 0.692, 0.242, 0.982, 0.264, 0.781, 0.779, 0.666, 0.549, 0.646, 0.05, 0.853, 0.729, 0.128, 0.111, 0.399, 0.353, 0.923, 0.962, 0.677, 0.324, 0.54, 0.678, 0.502, 0.475, 0.397, 0.942, 0.248, 0.963, 0.918, 0.769, 0.87, 0.436, 0.69, 0.737, 0.669, 0.839, 0.445, 0.961, 0.07, 0.885, 0.076, 0.311, 0.341, 0.048, 0.04, 0.65, 0.311, 0.437, 0.316, 0.652, 0.364, 0.875, 0.765, 0.758, 0.366, 0.274, 0.411, 0.918, 0.987, 0.06, 0.114, 0.277, 0.638, 0.778, 0.956, 0.086, 0.53, 0.902, 0.98, 0.517, 0.286, 0.497, 0.94, 0.119, 0.304, 0.129, 0.363, 0.529, 0.389, 0.498, 0.936, 0.248, 0.132, 0.109, 0.205, 0.243, 0.843, 0.674, 0.81, 0.857, 0.526, 0.115, 0.893, 0.92, 0.947, 0.808, 0.756, 0.509, 0.714, 0.606, 0.755, 0.927, 0.882, 0.6, 0.637, 0.318, 0.911, 0.475, 0.306, 0.148, 0.994, 0.586, 0.7, 0.704, 0.877, 0.227, 0.521, 0.388, 0.829, 0.739, 0.695, 0.575, 0.973, 0.043, 0.995, 0.535, 0.85, 0.786, 0.785, 0.099, 0.879, 0.653, 0.507, 0.39, 0.791, 0.941, 0.627, 0.509, 0.55, 0.448, 0.962, 0.048, 0.578, 0.047, 0.189, 0.287, 0.472, 0.582, 0.251, 0.964, 0.982, 0.998, 0.322, 0.926, 0.49, 0.995, 0.816, 0.928, 0.981, 0.924, 0.48, 0.97, 0.775, 0.424, 0.614, 0.354, 0.805, 0.584, 0.735, 0.649, 0.665, 0.568, 0.249, 0.065, 0.558, 0.952, 0.75]
global q = [0.73, 0.964, 0.875, 0.752, 0.705, 0.907, 0.895, 0.596, 0.292, 0.747, 0.925, 0.709, 0.754, 0.983, 0.749, 0.937, 0.643, 0.363, 0.819, 0.412, 0.831, 0.936, 0.514, 0.998, 0.452, 0.419, 0.953, 0.996, 0.987, 0.605, 0.13, 0.86, 0.94, 0.438, 0.972, 0.862, 0.804, 0.621, 0.726, 0.981, 0.912, 0.833, 0.802, 0.826, 0.272, 0.999, 0.999, 0.859, 0.539, 0.986, 0.632, 0.817, 0.931, 0.888, 0.901, 0.681, 0.06, 0.904, 0.93, 0.352, 0.697, 0.547, 0.888, 0.977, 0.978, 0.979, 0.873, 0.942, 0.932, 0.837, 0.658, 0.406, 0.997, 0.754, 0.969, 0.93, 0.984, 0.916, 0.452, 0.727, 0.858, 0.802, 0.88, 0.799, 0.99, 0.49, 0.889, 0.892, 0.338, 0.606, 0.234, 0.665, 0.995, 0.727, 0.974, 0.407, 0.72, 0.713, 0.994, 0.844, 0.774, 0.795, 0.624, 0.575, 0.925, 0.994, 0.607, 0.661, 0.613, 0.756, 0.909, 0.99, 0.705, 0.898, 0.985, 0.983, 0.98, 0.396, 0.706, 0.955, 0.784, 0.765, 0.628, 0.907, 0.691, 0.628, 0.564, 0.961, 0.253, 0.58, 0.438, 0.404, 0.687, 0.967, 0.814, 0.915, 0.955, 0.641, 0.229, 0.939, 0.988, 0.959, 0.947, 0.993, 0.816, 0.995, 0.972, 0.855, 0.938, 0.937, 0.875, 0.651, 0.879, 0.912, 0.952, 0.467, 0.599, 0.995, 0.643, 0.864, 0.742, 0.949, 0.747, 0.883, 0.641, 0.966, 0.932, 0.772, 0.686, 0.998, 0.53, 0.997, 0.913, 0.967, 0.88, 0.884, 0.241, 0.969, 0.784, 0.9, 0.715, 0.899, 0.966, 0.728, 0.917, 0.637, 0.577, 0.997, 0.611, 0.683, 0.652, 0.443, 0.969, 0.927, 0.83, 0.524, 0.994, 0.985, 0.998, 0.397, 0.937, 0.725, 0.997, 0.858, 0.982, 0.981, 0.98, 0.727, 0.97, 0.849, 0.673, 0.924, 0.764, 0.841, 0.842, 0.963, 0.895, 0.936, 0.707, 0.7, 0.106, 0.85, 0.956, 0.948]
global origin = 1
global destination = 50