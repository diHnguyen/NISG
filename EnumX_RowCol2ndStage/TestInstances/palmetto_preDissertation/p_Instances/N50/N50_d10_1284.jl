global arcs = [1 11; 1 27; 1 44; 2 30; 2 33; 2 47; 3 20; 3 22; 3 30; 3 35; 3 37; 4 7; 4 37; 4 47; 5 21; 5 31; 6 7; 6 24; 6 25; 6 45; 6 47; 7 2; 7 6; 7 18; 7 21; 7 28; 7 39; 8 3; 8 5; 8 31; 9 8; 9 10; 9 17; 9 26; 9 36; 9 48; 9 49; 9 50; 10 45; 11 9; 11 12; 11 18; 11 34; 11 41; 11 47; 12 24; 12 29; 12 46; 13 3; 13 32; 13 33; 14 4; 14 18; 14 38; 14 45; 15 14; 15 19; 15 22; 15 28; 15 35; 15 38; 15 40; 15 46; 16 4; 16 17; 16 28; 16 34; 16 39; 16 50; 17 8; 17 29; 17 41; 17 45; 18 2; 18 6; 18 9; 18 15; 18 17; 18 41; 18 49; 19 4; 19 5; 19 7; 19 9; 20 8; 20 22; 20 31; 20 41; 20 46; 21 17; 21 19; 21 22; 21 27; 21 35; 21 48; 22 6; 22 34; 23 5; 23 20; 23 22; 23 26; 23 37; 23 38; 23 40; 23 50; 24 5; 24 9; 24 13; 24 29; 24 37; 24 40; 24 43; 24 44; 25 3; 25 8; 25 20; 25 35; 25 36; 25 40; 26 16; 26 17; 26 25; 26 31; 26 38; 26 41; 27 6; 27 9; 27 20; 27 31; 28 3; 28 15; 28 21; 28 46; 29 3; 29 11; 29 19; 29 23; 29 24; 29 26; 30 7; 30 11; 30 12; 30 21; 30 29; 30 35; 30 48; 31 2; 31 6; 31 20; 31 23; 31 24; 31 29; 31 35; 31 44; 31 47; 32 3; 32 5; 32 9; 32 25; 32 43; 32 45; 33 10; 33 16; 33 21; 33 27; 33 38; 33 47; 34 2; 34 18; 34 42; 35 6; 35 10; 35 42; 36 18; 36 23; 36 24; 36 25; 36 37; 36 41; 36 43; 36 44; 37 10; 37 11; 37 45; 38 7; 38 14; 38 22; 38 33; 39 6; 39 12; 39 29; 39 34; 40 15; 40 17; 40 23; 40 24; 40 25; 40 37; 40 41; 41 12; 41 20; 41 21; 41 22; 41 28; 41 30; 41 33; 41 39; 41 50; 42 12; 42 15; 42 18; 42 31; 42 32; 42 37; 43 6; 43 31; 43 33; 43 36; 43 39; 44 26; 44 29; 44 32; 44 39; 45 3; 45 6; 45 33; 45 38; 46 2; 46 6; 46 14; 47 13; 47 20; 47 46; 48 3; 48 15; 48 29; 49 4; 49 8; 49 20; 49 37; 49 40]
global d_x = [4.0, 2.0, 1.0, 10.0, 5.0, 1.0, 3.0, 8.0, 5.0, 8.0, 6.0, 9.0, 5.0, 5.0, 9.0, 7.0, 2.0, 3.0, 2.0, 2.0, 9.0, 4.0, 4.0, 1.0, 2.0, 4.0, 10.0, 5.0, 5.0, 4.0, 10.0, 2.0, 1.0, 5.0, 3.0, 5.0, 10.0, 3.0, 2.0, 10.0, 1.0, 1.0, 1.0, 6.0, 2.0, 5.0, 3.0, 2.0, 7.0, 3.0, 7.0, 3.0, 2.0, 10.0, 8.0, 10.0, 4.0, 6.0, 9.0, 9.0, 4.0, 5.0, 3.0, 6.0, 4.0, 8.0, 4.0, 8.0, 6.0, 4.0, 2.0, 8.0, 8.0, 8.0, 3.0, 6.0, 7.0, 8.0, 7.0, 3.0, 5.0, 8.0, 10.0, 7.0, 10.0, 9.0, 7.0, 3.0, 10.0, 2.0, 7.0, 3.0, 7.0, 7.0, 7.0, 6.0, 10.0, 1.0, 9.0, 3.0, 7.0, 1.0, 2.0, 6.0, 4.0, 1.0, 2.0, 3.0, 6.0, 1.0, 4.0, 5.0, 10.0, 8.0, 1.0, 9.0, 5.0, 3.0, 7.0, 9.0, 1.0, 7.0, 1.0, 7.0, 4.0, 9.0, 2.0, 6.0, 10.0, 6.0, 1.0, 2.0, 5.0, 9.0, 6.0, 9.0, 2.0, 7.0, 6.0, 1.0, 3.0, 5.0, 7.0, 7.0, 4.0, 7.0, 7.0, 7.0, 8.0, 3.0, 4.0, 4.0, 5.0, 2.0, 8.0, 9.0, 5.0, 3.0, 8.0, 2.0, 8.0, 8.0, 5.0, 4.0, 6.0, 1.0, 4.0, 9.0, 4.0, 3.0, 4.0, 9.0, 6.0, 10.0, 4.0, 9.0, 1.0, 10.0, 3.0, 7.0, 3.0, 9.0, 3.0, 3.0, 5.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 3.0, 9.0, 6.0, 5.0, 1.0, 4.0, 3.0, 10.0, 8.0, 7.0, 7.0, 1.0, 4.0, 10.0, 10.0, 4.0, 6.0, 5.0, 1.0, 9.0, 8.0, 6.0, 2.0, 4.0, 2.0, 7.0, 6.0, 9.0, 1.0, 4.0, 5.0, 7.0, 3.0, 9.0, 6.0, 5.0, 1.0, 2.0, 6.0, 1.0, 2.0, 7.0, 5.0, 1.0, 6.0, 7.0, 4.0, 5.0, 7.0, 7.0]
global b_x = 5
global d_y = [7.0, 4.0, 7.0, 5.0, 8.0, 1.0, 6.0, 10.0, 5.0, 1.0, 5.0, 4.0, 2.0, 8.0, 8.0, 4.0, 9.0, 7.0, 6.0, 2.0, 8.0, 4.0, 5.0, 5.0, 4.0, 1.0, 10.0, 3.0, 8.0, 10.0, 7.0, 6.0, 1.0, 5.0, 5.0, 4.0, 10.0, 2.0, 2.0, 10.0, 5.0, 5.0, 3.0, 2.0, 7.0, 4.0, 9.0, 5.0, 10.0, 1.0, 2.0, 8.0, 1.0, 4.0, 1.0, 10.0, 5.0, 7.0, 8.0, 1.0, 10.0, 6.0, 2.0, 9.0, 3.0, 7.0, 9.0, 2.0, 6.0, 5.0, 9.0, 8.0, 2.0, 1.0, 2.0, 3.0, 9.0, 4.0, 2.0, 9.0, 3.0, 8.0, 9.0, 5.0, 1.0, 6.0, 4.0, 10.0, 6.0, 3.0, 7.0, 5.0, 4.0, 1.0, 7.0, 7.0, 4.0, 1.0, 1.0, 7.0, 5.0, 2.0, 8.0, 6.0, 7.0, 10.0, 9.0, 7.0, 1.0, 4.0, 7.0, 5.0, 10.0, 2.0, 9.0, 10.0, 8.0, 6.0, 7.0, 5.0, 7.0, 5.0, 8.0, 5.0, 1.0, 7.0, 9.0, 3.0, 7.0, 9.0, 9.0, 6.0, 3.0, 1.0, 10.0, 7.0, 9.0, 4.0, 9.0, 6.0, 10.0, 4.0, 10.0, 7.0, 9.0, 10.0, 3.0, 6.0, 7.0, 6.0, 9.0, 10.0, 9.0, 7.0, 9.0, 7.0, 8.0, 1.0, 9.0, 8.0, 6.0, 10.0, 9.0, 5.0, 3.0, 9.0, 1.0, 9.0, 7.0, 9.0, 9.0, 2.0, 6.0, 1.0, 3.0, 4.0, 5.0, 9.0, 3.0, 5.0, 6.0, 2.0, 4.0, 6.0, 4.0, 6.0, 3.0, 3.0, 4.0, 2.0, 2.0, 6.0, 3.0, 3.0, 9.0, 3.0, 10.0, 6.0, 6.0, 6.0, 1.0, 3.0, 6.0, 9.0, 5.0, 1.0, 2.0, 1.0, 1.0, 10.0, 6.0, 5.0, 8.0, 3.0, 2.0, 5.0, 9.0, 9.0, 2.0, 1.0, 2.0, 4.0, 2.0, 2.0, 9.0, 10.0, 9.0, 2.0, 9.0, 8.0, 10.0, 7.0, 7.0, 6.0, 3.0, 1.0, 2.0, 4.0, 7.0, 10.0, 5.0]
global b_y = 10
global p = [0.061, 0.73, 0.528, 0.997, 0.945, 0.272, 0.538, 0.064, 0.974, 0.183, 0.062, 0.333, 0.174, 0.211, 0.24, 0.33, 0.717, 0.808, 0.873, 0.028, 0.042, 0.796, 0.06, 0.782, 0.382, 0.428, 0.285, 0.192, 0.317, 0.201, 0.138, 0.522, 0.435, 0.345, 0.93, 0.356, 0.66, 0.872, 0.022, 0.249, 0.636, 0.449, 0.826, 0.672, 0.123, 0.856, 0.669, 0.666, 0.936, 0.908, 0.232, 0.861, 0.742, 0.107, 0.006, 0.251, 0.8, 0.61, 0.323, 0.152, 0.79, 0.923, 0.982, 0.747, 0.393, 0.674, 0.585, 0.378, 0.152, 0.847, 0.881, 0.695, 0.967, 0.745, 0.743, 0.506, 0.183, 0.14, 0.959, 0.214, 0.375, 0.597, 0.773, 0.034, 0.577, 0.677, 0.531, 0.951, 0.505, 0.318, 0.804, 0.367, 0.989, 0.262, 0.58, 0.088, 0.81, 0.756, 0.716, 0.225, 0.723, 0.536, 0.214, 0.142, 0.978, 0.051, 0.331, 0.474, 0.646, 0.548, 0.826, 0.707, 0.039, 0.094, 0.853, 0.874, 0.659, 0.03, 0.909, 0.339, 0.426, 0.21, 0.991, 0.549, 0.112, 0.435, 0.412, 0.202, 0.128, 0.026, 0.37, 0.544, 0.963, 0.463, 0.69, 0.593, 0.91, 0.781, 0.405, 0.445, 0.985, 0.121, 0.297, 0.579, 0.34, 0.837, 0.603, 0.17, 0.417, 0.812, 0.835, 0.194, 0.23, 0.926, 0.901, 0.205, 0.346, 0.677, 0.98, 0.643, 0.192, 0.893, 0.376, 0.299, 0.565, 0.579, 0.698, 0.041, 0.866, 0.515, 0.905, 0.178, 0.205, 0.172, 0.109, 0.041, 0.03, 0.637, 0.481, 0.281, 0.936, 0.24, 0.546, 0.351, 0.006, 0.301, 0.756, 0.34, 0.325, 0.361, 0.112, 0.026, 0.224, 0.071, 0.075, 0.547, 0.345, 0.165, 0.498, 0.727, 0.815, 0.046, 0.227, 0.618, 0.58, 0.953, 0.767, 0.514, 0.812, 0.029, 0.226, 0.858, 0.999, 0.284, 0.781, 0.027, 0.016, 0.803, 0.813, 0.87, 0.46, 0.007, 0.413, 0.21, 0.428, 0.204, 0.074, 0.176, 0.155, 0.427, 0.577, 0.928, 0.154, 0.64, 0.105, 0.047, 0.352, 0.674, 0.785, 0.63, 0.055]
global q = [0.838, 0.796, 0.56, 0.997, 0.984, 0.671, 0.595, 0.547, 0.977, 0.532, 0.777, 0.921, 0.725, 0.718, 0.992, 0.673, 0.843, 0.95, 0.905, 0.32, 0.862, 0.9, 0.862, 0.94, 0.736, 0.58, 0.965, 0.626, 0.596, 0.701, 0.588, 0.862, 0.455, 0.923, 0.971, 0.664, 0.988, 0.979, 0.171, 0.901, 0.69, 0.987, 0.834, 0.891, 0.168, 0.964, 0.705, 0.899, 0.937, 0.984, 0.55, 0.913, 0.915, 0.6, 0.359, 0.885, 0.882, 0.81, 0.449, 0.64, 0.859, 0.981, 0.982, 0.782, 0.531, 0.72, 0.783, 0.547, 0.304, 0.918, 0.945, 0.947, 0.984, 0.975, 0.872, 0.724, 0.377, 0.85, 0.986, 0.504, 0.881, 0.817, 0.962, 0.923, 0.882, 0.879, 0.659, 0.966, 0.813, 0.759, 0.912, 0.857, 0.991, 0.813, 0.916, 0.249, 0.829, 0.776, 0.781, 0.925, 0.835, 0.862, 0.357, 0.533, 0.988, 0.491, 0.442, 0.964, 0.752, 0.787, 0.839, 0.867, 0.425, 0.323, 0.966, 0.916, 0.803, 0.083, 0.92, 0.388, 0.466, 0.547, 0.995, 0.69, 0.494, 0.899, 0.599, 0.404, 0.482, 0.422, 0.811, 0.911, 0.964, 0.637, 0.942, 0.646, 0.929, 0.913, 0.778, 0.711, 0.998, 0.413, 0.428, 0.798, 0.421, 0.994, 0.659, 0.742, 0.936, 0.905, 0.864, 0.964, 0.885, 0.965, 0.96, 0.776, 0.922, 0.786, 0.986, 0.787, 0.771, 0.916, 0.506, 0.84, 0.915, 0.734, 0.838, 0.896, 0.998, 0.858, 0.938, 0.361, 0.652, 0.805, 0.798, 0.862, 0.952, 0.665, 0.78, 0.954, 0.952, 0.29, 0.823, 0.547, 0.132, 0.461, 0.996, 0.997, 0.47, 0.419, 0.779, 0.84, 0.313, 0.112, 0.894, 0.724, 0.765, 0.412, 0.593, 0.968, 0.897, 0.129, 0.333, 0.975, 0.775, 0.967, 0.883, 0.687, 0.836, 0.492, 0.897, 0.906, 0.999, 0.39, 0.886, 0.797, 0.343, 0.886, 0.956, 0.928, 0.769, 0.46, 0.791, 0.249, 0.83, 0.372, 0.796, 0.513, 0.642, 0.726, 0.874, 0.96, 0.761, 0.944, 0.702, 0.18, 0.781, 0.857, 0.908, 0.759, 0.907]
global origin = 1
global destination = 50