global arcs = [1 7; 1 8; 1 19; 1 25; 1 26; 2 15; 2 24; 3 8; 3 15; 3 20; 3 29; 3 30; 3 37; 4 29; 4 32; 4 36; 5 11; 5 15; 5 16; 5 19; 5 22; 6 4; 6 20; 6 27; 7 14; 7 20; 7 29; 8 9; 8 10; 9 12; 9 31; 9 36; 10 11; 10 17; 10 32; 11 9; 11 10; 11 21; 11 26; 11 29; 11 39; 11 40; 12 7; 12 22; 12 23; 12 36; 13 3; 13 10; 13 14; 13 20; 13 27; 13 28; 13 34; 14 7; 14 28; 14 32; 14 40; 15 2; 15 12; 15 13; 15 22; 15 25; 15 27; 15 33; 15 34; 15 38; 15 40; 16 6; 16 23; 16 26; 17 15; 17 25; 17 27; 17 29; 18 10; 18 22; 18 31; 18 32; 18 34; 19 17; 19 20; 19 22; 19 33; 19 36; 19 37; 20 9; 20 12; 20 19; 20 24; 20 25; 20 27; 20 36; 21 11; 21 15; 21 16; 21 35; 21 37; 22 33; 22 39; 23 24; 23 35; 24 5; 24 19; 24 27; 24 31; 24 34; 24 38; 25 17; 25 21; 25 33; 25 35; 26 2; 26 13; 26 20; 26 29; 26 33; 26 40; 27 4; 27 16; 27 32; 28 4; 28 5; 28 7; 28 12; 28 23; 29 10; 29 26; 29 28; 30 33; 31 26; 31 30; 32 14; 32 16; 32 18; 32 28; 32 35; 32 37; 33 5; 33 12; 33 20; 34 13; 34 17; 34 22; 35 17; 35 18; 35 21; 35 31; 35 38; 36 3; 36 8; 36 12; 36 13; 37 4; 37 6; 37 12; 37 20; 37 36; 38 2; 38 5; 38 8; 38 14; 38 37; 38 39; 39 3; 39 5; 39 23; 39 33]
global d_x = [3.0, 10.0, 3.0, 1.0, 4.0, 2.0, 3.0, 5.0, 4.0, 4.0, 6.0, 10.0, 6.0, 4.0, 8.0, 5.0, 1.0, 6.0, 6.0, 4.0, 7.0, 3.0, 4.0, 1.0, 6.0, 2.0, 4.0, 1.0, 9.0, 2.0, 8.0, 9.0, 5.0, 7.0, 8.0, 1.0, 7.0, 3.0, 1.0, 5.0, 10.0, 2.0, 2.0, 10.0, 4.0, 3.0, 7.0, 7.0, 2.0, 6.0, 2.0, 2.0, 8.0, 2.0, 10.0, 3.0, 10.0, 4.0, 8.0, 8.0, 1.0, 1.0, 2.0, 4.0, 8.0, 10.0, 9.0, 8.0, 4.0, 1.0, 3.0, 5.0, 4.0, 9.0, 5.0, 8.0, 9.0, 7.0, 9.0, 3.0, 1.0, 3.0, 7.0, 2.0, 7.0, 3.0, 3.0, 7.0, 3.0, 8.0, 6.0, 6.0, 6.0, 10.0, 9.0, 7.0, 9.0, 7.0, 5.0, 5.0, 1.0, 3.0, 6.0, 3.0, 3.0, 10.0, 1.0, 1.0, 10.0, 1.0, 4.0, 1.0, 10.0, 2.0, 9.0, 1.0, 8.0, 5.0, 5.0, 8.0, 4.0, 7.0, 7.0, 1.0, 8.0, 3.0, 9.0, 6.0, 1.0, 4.0, 5.0, 3.0, 1.0, 9.0, 4.0, 10.0, 8.0, 4.0, 2.0, 8.0, 9.0, 10.0, 10.0, 9.0, 9.0, 10.0, 7.0, 7.0, 3.0, 6.0, 4.0, 7.0, 3.0, 7.0, 7.0, 10.0, 5.0, 8.0, 6.0, 1.0, 9.0, 2.0, 4.0, 7.0, 3.0, 7.0, 10.0]
global b_x = 5
global d_y = [9.0, 5.0, 3.0, 5.0, 7.0, 2.0, 4.0, 6.0, 9.0, 3.0, 9.0, 2.0, 10.0, 5.0, 2.0, 1.0, 7.0, 10.0, 1.0, 6.0, 9.0, 1.0, 2.0, 9.0, 9.0, 7.0, 7.0, 1.0, 10.0, 4.0, 7.0, 7.0, 4.0, 9.0, 6.0, 7.0, 7.0, 4.0, 3.0, 6.0, 9.0, 7.0, 8.0, 10.0, 4.0, 5.0, 1.0, 5.0, 1.0, 8.0, 2.0, 6.0, 8.0, 9.0, 3.0, 9.0, 3.0, 6.0, 5.0, 8.0, 4.0, 3.0, 4.0, 3.0, 9.0, 9.0, 3.0, 4.0, 6.0, 4.0, 5.0, 1.0, 7.0, 1.0, 1.0, 5.0, 10.0, 2.0, 6.0, 4.0, 1.0, 2.0, 9.0, 2.0, 1.0, 5.0, 5.0, 5.0, 9.0, 8.0, 10.0, 5.0, 7.0, 2.0, 6.0, 1.0, 9.0, 7.0, 10.0, 8.0, 2.0, 2.0, 1.0, 7.0, 2.0, 2.0, 10.0, 9.0, 10.0, 3.0, 10.0, 3.0, 7.0, 6.0, 5.0, 2.0, 1.0, 8.0, 2.0, 7.0, 1.0, 8.0, 6.0, 9.0, 6.0, 1.0, 6.0, 3.0, 4.0, 5.0, 3.0, 2.0, 4.0, 4.0, 7.0, 1.0, 4.0, 2.0, 10.0, 7.0, 3.0, 5.0, 10.0, 8.0, 3.0, 10.0, 5.0, 1.0, 8.0, 1.0, 4.0, 3.0, 8.0, 5.0, 6.0, 6.0, 10.0, 2.0, 8.0, 1.0, 8.0, 2.0, 10.0, 10.0, 10.0, 9.0, 4.0]
global b_y = 10
global p = [0.586, 0.238, 0.258, 0.923, 0.459, 0.009, 0.462, 0.666, 0.733, 0.326, 0.533, 0.358, 0.123, 0.355, 0.367, 0.427, 0.717, 0.398, 0.324, 0.527, 0.401, 0.488, 0.839, 0.689, 0.102, 0.227, 0.442, 0.054, 0.603, 0.631, 0.676, 0.711, 0.745, 0.376, 0.562, 0.652, 0.25, 0.64, 0.526, 0.46, 0.128, 0.597, 0.059, 0.888, 0.538, 0.856, 0.678, 0.512, 0.218, 0.132, 0.551, 0.862, 0.344, 0.86, 0.102, 0.846, 0.259, 0.426, 0.864, 0.346, 0.85, 0.427, 0.267, 0.352, 0.326, 0.689, 0.708, 0.888, 0.88, 0.782, 0.8, 0.514, 0.581, 0.019, 0.123, 0.972, 0.201, 0.493, 0.045, 0.19, 0.068, 0.859, 0.958, 0.314, 0.798, 0.097, 0.332, 0.893, 0.683, 0.01, 0.516, 0.027, 0.327, 0.44, 0.087, 0.878, 0.454, 0.919, 0.295, 0.147, 0.304, 0.506, 0.395, 0.891, 0.358, 0.291, 0.375, 0.469, 0.841, 0.663, 0.419, 0.176, 0.595, 0.196, 0.66, 0.202, 0.656, 0.115, 0.137, 0.488, 0.483, 0.063, 0.538, 0.373, 0.917, 0.285, 0.071, 0.287, 0.763, 0.923, 0.074, 0.389, 0.079, 0.075, 0.208, 0.86, 0.02, 0.83, 0.162, 0.481, 0.47, 0.185, 0.133, 0.823, 0.924, 0.641, 0.686, 0.754, 0.437, 0.811, 0.59, 0.567, 0.129, 0.685, 0.626, 0.958, 0.84, 0.835, 0.691, 0.927, 0.525, 0.7, 0.937, 0.406, 0.897, 0.901, 0.465]
global q = [0.618, 0.263, 0.274, 0.986, 0.985, 0.492, 0.752, 0.871, 0.976, 0.861, 0.861, 0.814, 0.242, 0.624, 0.99, 0.434, 0.807, 0.562, 0.919, 0.865, 0.84, 0.987, 0.91, 0.801, 0.496, 0.689, 0.768, 0.175, 0.905, 0.869, 0.894, 0.939, 0.802, 0.435, 0.922, 0.829, 0.305, 0.865, 0.711, 0.495, 0.147, 0.807, 0.996, 0.957, 0.931, 0.896, 0.691, 0.78, 0.659, 0.976, 0.571, 0.939, 0.916, 0.887, 0.231, 0.941, 0.789, 0.463, 0.924, 0.776, 0.869, 0.491, 0.473, 0.556, 0.879, 0.997, 0.791, 0.963, 0.95, 0.865, 0.961, 0.623, 0.652, 0.778, 0.851, 0.982, 0.84, 0.72, 0.127, 0.872, 0.118, 0.928, 0.999, 0.693, 0.939, 0.755, 0.584, 0.911, 0.942, 0.188, 0.692, 0.974, 0.356, 0.495, 0.211, 0.888, 0.946, 0.987, 0.794, 0.694, 0.747, 0.938, 0.844, 0.969, 0.754, 0.361, 0.776, 0.696, 0.862, 0.861, 0.45, 0.442, 0.888, 0.627, 0.661, 0.333, 0.785, 0.646, 0.49, 0.635, 0.615, 0.602, 0.766, 0.864, 0.93, 0.835, 0.233, 0.877, 0.829, 0.959, 0.118, 0.612, 0.964, 0.995, 0.943, 0.939, 0.329, 0.899, 0.762, 0.516, 0.589, 0.248, 0.683, 0.998, 0.997, 0.811, 0.911, 0.962, 0.619, 0.895, 0.816, 0.944, 0.567, 0.98, 0.905, 0.994, 0.872, 0.934, 0.878, 0.944, 0.999, 0.943, 0.94, 0.467, 0.986, 0.986, 0.931]
global origin = 1
global destination = 40
