global arcs = [1 8; 1 11; 1 20; 1 28; 1 36; 1 40; 2 24; 2 26; 2 29; 2 32; 3 2; 3 22; 3 41; 3 44; 4 8; 4 9; 4 11; 4 23; 4 34; 4 37; 5 23; 6 2; 6 36; 6 44; 7 9; 7 11; 7 19; 7 26; 7 27; 7 41; 7 47; 8 23; 8 24; 8 35; 9 28; 9 30; 9 35; 9 41; 9 43; 10 17; 10 22; 10 35; 10 40; 10 42; 10 44; 10 45; 11 3; 11 32; 12 6; 12 18; 12 34; 13 23; 13 31; 13 32; 13 49; 14 3; 14 26; 14 50; 15 3; 15 11; 15 17; 15 36; 16 7; 16 13; 16 22; 17 2; 17 9; 17 13; 17 19; 17 38; 18 2; 18 10; 18 29; 18 43; 18 49; 19 2; 19 3; 19 5; 19 14; 19 26; 19 36; 19 42; 19 45; 19 46; 19 47; 19 50; 20 2; 20 3; 20 5; 20 9; 20 28; 20 43; 20 44; 21 6; 21 48; 22 9; 22 13; 22 25; 22 32; 22 37; 22 48; 23 6; 23 7; 23 9; 23 19; 23 22; 23 24; 23 30; 23 33; 23 34; 23 37; 23 40; 24 3; 24 15; 24 29; 25 10; 25 24; 25 33; 25 34; 25 41; 25 43; 25 46; 26 4; 26 5; 26 16; 26 24; 26 34; 26 42; 26 43; 26 44; 27 10; 27 25; 28 2; 28 24; 28 25; 28 33; 29 18; 29 28; 29 33; 29 34; 29 45; 30 8; 30 11; 30 24; 30 27; 30 35; 30 43; 31 18; 31 21; 31 23; 32 3; 32 14; 32 20; 32 21; 32 27; 33 7; 33 13; 33 20; 33 40; 33 45; 34 2; 34 11; 34 12; 35 11; 35 18; 35 22; 35 34; 35 36; 36 20; 36 41; 37 10; 37 11; 37 14; 37 17; 37 27; 37 43; 37 46; 37 50; 38 3; 38 8; 38 10; 38 12; 38 25; 39 21; 39 24; 39 25; 39 31; 39 43; 40 9; 40 15; 40 25; 40 29; 41 2; 41 18; 41 20; 42 9; 42 34; 43 24; 43 33; 43 46; 43 48; 44 4; 44 27; 45 8; 45 22; 45 47; 45 49; 46 22; 46 39; 46 48; 47 3; 47 8; 47 19; 47 24; 47 36; 48 6; 48 9; 48 14; 48 16; 48 21; 49 5; 49 9; 49 19; 49 21; 49 29; 49 43]
global d_x = [4.0, 5.0, 8.0, 4.0, 10.0, 3.0, 9.0, 10.0, 8.0, 6.0, 6.0, 1.0, 3.0, 2.0, 6.0, 3.0, 10.0, 6.0, 1.0, 2.0, 9.0, 1.0, 9.0, 1.0, 2.0, 1.0, 9.0, 10.0, 4.0, 2.0, 5.0, 9.0, 7.0, 2.0, 6.0, 8.0, 9.0, 7.0, 1.0, 8.0, 6.0, 3.0, 6.0, 6.0, 2.0, 4.0, 6.0, 10.0, 3.0, 8.0, 9.0, 10.0, 3.0, 6.0, 1.0, 10.0, 7.0, 4.0, 5.0, 5.0, 1.0, 6.0, 1.0, 6.0, 1.0, 10.0, 3.0, 10.0, 2.0, 3.0, 4.0, 9.0, 8.0, 3.0, 9.0, 2.0, 10.0, 4.0, 6.0, 8.0, 1.0, 5.0, 10.0, 8.0, 2.0, 7.0, 9.0, 6.0, 4.0, 7.0, 5.0, 9.0, 3.0, 9.0, 9.0, 7.0, 2.0, 5.0, 2.0, 2.0, 2.0, 6.0, 2.0, 7.0, 5.0, 1.0, 2.0, 1.0, 10.0, 9.0, 2.0, 2.0, 2.0, 5.0, 4.0, 8.0, 7.0, 3.0, 1.0, 7.0, 2.0, 8.0, 6.0, 8.0, 7.0, 6.0, 8.0, 3.0, 2.0, 4.0, 1.0, 1.0, 9.0, 7.0, 3.0, 10.0, 1.0, 10.0, 5.0, 2.0, 8.0, 3.0, 6.0, 3.0, 7.0, 10.0, 5.0, 7.0, 1.0, 9.0, 5.0, 7.0, 3.0, 8.0, 5.0, 1.0, 5.0, 8.0, 2.0, 6.0, 1.0, 1.0, 6.0, 9.0, 6.0, 9.0, 1.0, 9.0, 9.0, 10.0, 7.0, 2.0, 3.0, 7.0, 7.0, 1.0, 4.0, 7.0, 6.0, 10.0, 8.0, 2.0, 9.0, 1.0, 1.0, 6.0, 2.0, 7.0, 2.0, 4.0, 6.0, 3.0, 10.0, 2.0, 2.0, 3.0, 7.0, 1.0, 4.0, 10.0, 3.0, 4.0, 8.0, 2.0, 2.0, 9.0, 5.0, 4.0, 7.0, 1.0, 3.0, 9.0, 4.0, 10.0, 1.0, 4.0, 7.0, 10.0, 8.0, 10.0, 8.0, 2.0, 4.0, 9.0, 7.0, 9.0]
global b_x = 5
global d_y = [8.0, 4.0, 8.0, 6.0, 1.0, 7.0, 1.0, 3.0, 5.0, 2.0, 9.0, 9.0, 10.0, 7.0, 4.0, 1.0, 8.0, 6.0, 4.0, 5.0, 3.0, 8.0, 6.0, 1.0, 9.0, 6.0, 1.0, 7.0, 7.0, 3.0, 6.0, 3.0, 7.0, 9.0, 7.0, 5.0, 3.0, 4.0, 7.0, 6.0, 5.0, 1.0, 2.0, 1.0, 7.0, 5.0, 6.0, 4.0, 4.0, 4.0, 7.0, 10.0, 3.0, 8.0, 8.0, 2.0, 3.0, 9.0, 3.0, 1.0, 2.0, 4.0, 1.0, 5.0, 2.0, 10.0, 2.0, 5.0, 7.0, 4.0, 1.0, 10.0, 3.0, 8.0, 3.0, 10.0, 9.0, 3.0, 10.0, 10.0, 4.0, 1.0, 5.0, 10.0, 10.0, 8.0, 8.0, 3.0, 3.0, 1.0, 1.0, 8.0, 4.0, 1.0, 4.0, 2.0, 1.0, 6.0, 6.0, 7.0, 10.0, 3.0, 8.0, 2.0, 10.0, 2.0, 6.0, 2.0, 1.0, 10.0, 6.0, 8.0, 5.0, 4.0, 3.0, 8.0, 6.0, 2.0, 7.0, 3.0, 4.0, 7.0, 10.0, 2.0, 8.0, 10.0, 1.0, 9.0, 10.0, 10.0, 3.0, 5.0, 10.0, 2.0, 1.0, 7.0, 9.0, 5.0, 2.0, 6.0, 9.0, 6.0, 4.0, 3.0, 3.0, 9.0, 8.0, 1.0, 7.0, 1.0, 6.0, 5.0, 9.0, 7.0, 6.0, 7.0, 1.0, 4.0, 3.0, 2.0, 9.0, 3.0, 7.0, 7.0, 2.0, 7.0, 2.0, 6.0, 7.0, 1.0, 7.0, 7.0, 4.0, 7.0, 9.0, 9.0, 8.0, 7.0, 5.0, 9.0, 3.0, 5.0, 7.0, 7.0, 8.0, 4.0, 7.0, 4.0, 4.0, 9.0, 7.0, 9.0, 1.0, 4.0, 1.0, 3.0, 10.0, 9.0, 10.0, 3.0, 4.0, 8.0, 3.0, 10.0, 10.0, 10.0, 3.0, 6.0, 8.0, 4.0, 7.0, 2.0, 6.0, 1.0, 3.0, 10.0, 6.0, 3.0, 6.0, 4.0, 2.0, 10.0, 2.0, 3.0, 8.0, 9.0]
global b_y = 10
global p = [0.204, 0.682, 0.319, 0.184, 0.463, 0.756, 0.724, 0.226, 0.898, 0.599, 0.211, 0.609, 0.057, 0.41, 0.728, 0.455, 0.287, 0.719, 0.41, 0.948, 0.916, 0.931, 0.341, 0.373, 0.752, 0.311, 0.893, 0.64, 0.301, 0.949, 0.621, 0.239, 0.925, 0.669, 0.062, 0.6, 0.316, 0.298, 0.882, 0.079, 0.531, 0.423, 0.442, 0.484, 0.915, 0.808, 0.057, 0.729, 0.835, 0.968, 0.551, 0.255, 0.711, 0.318, 0.391, 0.276, 0.807, 0.016, 0.267, 0.994, 0.26, 0.813, 0.104, 0.562, 0.772, 0.862, 0.944, 0.022, 0.232, 0.034, 0.134, 0.116, 0.077, 0.354, 0.591, 0.658, 0.48, 0.493, 0.19, 0.455, 0.411, 0.228, 0.784, 0.551, 0.055, 0.636, 0.874, 0.952, 0.256, 0.048, 0.2, 0.91, 0.006, 0.217, 0.371, 0.327, 0.848, 0.671, 0.287, 0.213, 0.294, 0.586, 0.908, 0.956, 0.462, 0.724, 0.34, 0.541, 0.887, 0.536, 0.148, 0.405, 0.348, 0.377, 0.537, 0.335, 0.436, 0.336, 0.736, 0.551, 0.262, 0.954, 0.208, 0.886, 0.45, 0.593, 0.614, 0.909, 0.551, 0.402, 0.387, 0.917, 0.844, 0.666, 0.29, 0.551, 0.647, 0.078, 0.354, 0.609, 0.868, 0.355, 0.406, 0.974, 0.219, 0.254, 0.266, 0.245, 0.426, 0.899, 0.036, 0.163, 0.508, 0.39, 0.107, 0.257, 0.089, 0.713, 0.768, 0.994, 0.214, 0.34, 0.766, 0.105, 0.111, 0.509, 0.425, 0.023, 0.316, 0.128, 0.31, 0.42, 0.392, 0.495, 0.643, 0.471, 0.953, 0.485, 0.941, 0.909, 0.204, 0.83, 0.883, 0.551, 0.61, 0.954, 0.634, 0.942, 0.704, 0.922, 0.266, 0.074, 0.139, 0.829, 0.47, 0.298, 0.529, 0.245, 0.974, 0.34, 0.817, 0.182, 0.433, 0.98, 0.613, 0.843, 0.585, 0.786, 0.241, 0.434, 0.417, 0.597, 0.65, 0.654, 0.592, 0.158, 0.662, 0.018, 0.735, 0.189, 0.728, 0.118, 0.398, 0.264, 0.439, 0.755]
global q = [0.843, 0.825, 0.65, 0.58, 0.832, 0.882, 0.966, 0.874, 0.898, 0.853, 0.483, 0.93, 0.343, 0.446, 0.778, 0.972, 0.957, 0.812, 0.722, 0.959, 0.998, 0.975, 0.471, 0.639, 0.764, 0.747, 0.922, 0.876, 0.735, 0.987, 0.921, 0.551, 0.992, 0.829, 0.53, 0.963, 0.975, 0.899, 0.902, 0.51, 0.607, 0.619, 0.751, 0.8, 0.928, 0.904, 0.667, 0.742, 0.935, 0.99, 0.591, 0.923, 0.829, 0.401, 0.51, 0.933, 0.969, 0.581, 0.523, 0.996, 0.602, 0.906, 0.996, 0.955, 0.996, 0.999, 0.962, 0.231, 0.895, 0.767, 0.204, 0.628, 0.826, 0.706, 0.817, 0.659, 0.511, 0.554, 0.75, 0.938, 0.442, 0.999, 0.897, 0.597, 0.967, 0.795, 0.961, 0.965, 0.441, 0.701, 0.424, 0.923, 0.937, 0.391, 0.971, 0.78, 0.992, 0.844, 0.477, 0.292, 0.791, 0.9, 0.918, 0.995, 0.911, 0.725, 0.446, 0.86, 0.936, 0.884, 0.819, 0.762, 0.943, 0.809, 0.638, 0.677, 0.759, 0.55, 0.993, 0.555, 0.55, 0.985, 0.351, 0.998, 0.693, 0.743, 0.969, 0.969, 0.849, 0.783, 0.74, 0.941, 0.989, 0.914, 0.599, 0.98, 0.684, 0.949, 0.727, 0.786, 0.997, 0.477, 0.996, 0.989, 0.354, 0.87, 0.893, 0.865, 0.528, 0.99, 0.547, 0.587, 0.967, 0.4, 0.285, 0.755, 0.629, 0.82, 0.854, 0.995, 0.831, 0.941, 0.781, 0.948, 0.529, 0.842, 0.522, 0.483, 0.458, 0.499, 0.659, 0.968, 0.853, 0.858, 0.685, 0.573, 0.977, 0.891, 0.995, 0.927, 0.332, 0.889, 0.953, 0.644, 0.956, 0.975, 0.997, 0.948, 0.978, 0.951, 0.548, 0.201, 0.514, 0.897, 0.781, 0.78, 0.627, 0.651, 0.997, 0.794, 0.835, 0.931, 0.695, 0.997, 0.646, 0.908, 0.846, 0.813, 0.787, 0.718, 0.612, 0.9, 0.744, 0.685, 0.607, 0.78, 0.843, 0.775, 0.736, 0.803, 0.902, 0.822, 0.836, 0.948, 0.976, 0.765]
global origin = 1
global destination = 50