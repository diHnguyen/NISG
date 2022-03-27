global arcs = [1 6; 1 16; 1 17; 1 36; 1 40; 1 43; 2 5; 2 31; 2 33; 2 39; 2 46; 2 47; 3 20; 3 24; 3 32; 3 36; 3 41; 4 13; 4 14; 4 20; 4 31; 4 32; 4 41; 4 43; 5 6; 5 24; 5 33; 5 48; 6 7; 6 26; 6 27; 6 34; 6 46; 6 48; 7 23; 7 26; 8 17; 8 18; 8 33; 8 34; 8 37; 8 48; 9 6; 9 27; 10 13; 10 16; 10 23; 10 42; 11 15; 11 36; 11 39; 11 48; 12 18; 12 24; 12 38; 12 48; 12 50; 13 24; 13 31; 14 19; 14 23; 14 30; 14 33; 14 44; 14 46; 14 50; 15 20; 15 42; 16 21; 16 41; 16 47; 17 7; 17 45; 18 5; 18 10; 18 21; 18 47; 19 3; 19 29; 19 30; 19 40; 20 2; 20 34; 20 39; 21 10; 21 25; 21 35; 21 36; 21 40; 21 50; 22 3; 22 7; 22 10; 22 11; 22 20; 22 27; 22 47; 22 48; 23 11; 23 15; 23 18; 23 25; 24 2; 24 7; 24 19; 24 32; 25 5; 25 9; 25 10; 25 24; 25 27; 25 40; 25 44; 26 8; 26 10; 26 29; 26 40; 26 49; 27 3; 27 5; 27 12; 27 17; 27 36; 28 13; 28 18; 28 24; 28 39; 29 4; 29 19; 29 21; 29 24; 29 25; 29 36; 29 43; 30 5; 30 19; 30 28; 30 32; 30 38; 31 5; 31 27; 31 47; 32 2; 32 11; 32 15; 32 33; 32 35; 32 48; 33 14; 33 21; 33 44; 34 11; 34 41; 35 4; 35 5; 35 7; 35 19; 35 30; 35 50; 36 3; 36 5; 36 7; 36 17; 36 18; 36 26; 36 32; 37 11; 37 44; 37 45; 37 47; 37 48; 38 35; 38 43; 38 48; 39 2; 39 14; 39 31; 39 43; 39 45; 39 50; 40 9; 40 25; 40 28; 40 30; 40 33; 40 35; 40 38; 40 48; 41 11; 41 21; 41 26; 41 43; 42 4; 42 6; 42 14; 42 22; 42 33; 42 38; 42 45; 42 48; 42 49; 42 50; 43 14; 43 27; 44 7; 44 21; 44 34; 45 21; 45 46; 46 2; 46 3; 46 5; 46 16; 46 17; 46 20; 46 29; 46 33; 46 45; 46 49; 47 22; 47 45; 48 8; 48 29; 48 30; 48 34; 49 3; 49 21; 49 31]
global d_x = [7.0, 2.0, 9.0, 7.0, 7.0, 8.0, 3.0, 8.0, 7.0, 5.0, 5.0, 1.0, 5.0, 3.0, 1.0, 9.0, 1.0, 2.0, 3.0, 1.0, 8.0, 2.0, 4.0, 8.0, 6.0, 5.0, 10.0, 4.0, 2.0, 8.0, 10.0, 9.0, 5.0, 2.0, 5.0, 6.0, 8.0, 3.0, 1.0, 3.0, 6.0, 3.0, 7.0, 7.0, 3.0, 5.0, 2.0, 8.0, 9.0, 3.0, 5.0, 2.0, 4.0, 5.0, 2.0, 1.0, 3.0, 5.0, 4.0, 4.0, 5.0, 4.0, 10.0, 8.0, 5.0, 9.0, 5.0, 5.0, 4.0, 3.0, 10.0, 8.0, 5.0, 8.0, 6.0, 1.0, 10.0, 6.0, 8.0, 4.0, 6.0, 10.0, 8.0, 10.0, 9.0, 1.0, 7.0, 7.0, 10.0, 5.0, 9.0, 10.0, 5.0, 2.0, 7.0, 10.0, 5.0, 2.0, 10.0, 4.0, 1.0, 3.0, 9.0, 3.0, 4.0, 6.0, 2.0, 5.0, 7.0, 3.0, 1.0, 10.0, 9.0, 3.0, 9.0, 2.0, 4.0, 3.0, 3.0, 1.0, 6.0, 2.0, 2.0, 7.0, 4.0, 3.0, 2.0, 7.0, 2.0, 10.0, 3.0, 10.0, 5.0, 1.0, 7.0, 5.0, 3.0, 6.0, 10.0, 10.0, 2.0, 9.0, 3.0, 1.0, 8.0, 6.0, 1.0, 10.0, 2.0, 8.0, 7.0, 5.0, 7.0, 7.0, 9.0, 4.0, 2.0, 3.0, 8.0, 2.0, 2.0, 9.0, 3.0, 1.0, 9.0, 8.0, 4.0, 8.0, 6.0, 2.0, 1.0, 5.0, 2.0, 8.0, 4.0, 5.0, 1.0, 3.0, 10.0, 4.0, 2.0, 10.0, 7.0, 9.0, 7.0, 6.0, 4.0, 4.0, 2.0, 9.0, 3.0, 5.0, 10.0, 10.0, 7.0, 2.0, 9.0, 7.0, 4.0, 5.0, 3.0, 7.0, 2.0, 10.0, 7.0, 3.0, 3.0, 4.0, 3.0, 10.0, 10.0, 9.0, 5.0, 10.0, 3.0, 7.0, 6.0, 8.0, 7.0, 4.0, 3.0, 1.0, 4.0, 9.0, 4.0, 8.0, 10.0, 4.0]
global b_x = 5
global d_y = [7.0, 10.0, 6.0, 4.0, 5.0, 8.0, 6.0, 4.0, 2.0, 10.0, 10.0, 8.0, 1.0, 5.0, 8.0, 6.0, 1.0, 8.0, 8.0, 5.0, 8.0, 3.0, 5.0, 5.0, 7.0, 5.0, 1.0, 6.0, 1.0, 9.0, 1.0, 10.0, 3.0, 7.0, 8.0, 8.0, 8.0, 5.0, 3.0, 9.0, 1.0, 8.0, 10.0, 3.0, 4.0, 7.0, 5.0, 5.0, 7.0, 7.0, 6.0, 8.0, 4.0, 2.0, 6.0, 2.0, 10.0, 1.0, 5.0, 6.0, 6.0, 6.0, 1.0, 6.0, 3.0, 9.0, 1.0, 1.0, 7.0, 4.0, 2.0, 4.0, 9.0, 10.0, 4.0, 5.0, 8.0, 10.0, 3.0, 6.0, 1.0, 5.0, 2.0, 8.0, 10.0, 8.0, 2.0, 7.0, 2.0, 7.0, 7.0, 3.0, 5.0, 4.0, 10.0, 4.0, 7.0, 1.0, 3.0, 1.0, 8.0, 4.0, 6.0, 6.0, 4.0, 1.0, 5.0, 9.0, 7.0, 6.0, 1.0, 1.0, 1.0, 3.0, 7.0, 5.0, 8.0, 10.0, 5.0, 10.0, 9.0, 6.0, 3.0, 6.0, 6.0, 9.0, 8.0, 2.0, 9.0, 9.0, 4.0, 4.0, 5.0, 5.0, 5.0, 6.0, 9.0, 4.0, 7.0, 4.0, 7.0, 8.0, 9.0, 9.0, 2.0, 6.0, 7.0, 10.0, 10.0, 9.0, 4.0, 7.0, 3.0, 10.0, 2.0, 2.0, 1.0, 2.0, 10.0, 2.0, 8.0, 8.0, 5.0, 3.0, 2.0, 1.0, 8.0, 4.0, 10.0, 3.0, 8.0, 5.0, 2.0, 1.0, 9.0, 1.0, 5.0, 9.0, 9.0, 5.0, 2.0, 6.0, 8.0, 10.0, 5.0, 6.0, 3.0, 6.0, 3.0, 7.0, 1.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 1.0, 9.0, 2.0, 8.0, 6.0, 7.0, 6.0, 10.0, 5.0, 8.0, 9.0, 2.0, 1.0, 3.0, 3.0, 9.0, 6.0, 7.0, 2.0, 5.0, 5.0, 2.0, 1.0, 6.0, 4.0, 5.0, 7.0, 8.0, 3.0, 6.0, 7.0]
global b_y = 10
global p = [0.028, 0.999, 0.979, 0.202, 0.998, 0.943, 0.888, 0.604, 0.838, 0.479, 0.251, 0.729, 0.463, 0.578, 0.635, 0.689, 0.184, 0.117, 0.357, 0.609, 0.321, 0.477, 0.272, 0.363, 0.656, 0.71, 0.201, 0.209, 0.788, 0.665, 0.8, 0.203, 0.45, 0.413, 0.479, 0.864, 0.978, 0.351, 0.949, 0.003, 0.212, 0.885, 0.844, 0.514, 0.59, 0.199, 0.188, 0.658, 0.926, 0.931, 0.894, 0.471, 0.794, 0.038, 0.625, 0.321, 0.165, 0.644, 0.07, 0.778, 0.009, 0.05, 0.504, 0.994, 0.286, 0.572, 0.395, 0.502, 0.158, 0.155, 0.264, 0.302, 0.499, 0.707, 0.58, 0.649, 0.532, 0.638, 0.206, 0.539, 0.985, 0.08, 0.788, 0.79, 0.603, 0.747, 0.572, 0.158, 0.785, 0.054, 0.262, 0.5, 0.839, 0.728, 0.433, 0.464, 0.085, 0.953, 0.846, 0.99, 0.841, 0.458, 0.458, 0.626, 0.73, 0.961, 0.987, 0.248, 0.509, 0.341, 0.206, 0.973, 0.274, 0.059, 0.902, 0.785, 0.748, 0.029, 0.457, 0.68, 0.215, 0.103, 0.622, 0.377, 0.462, 0.023, 0.174, 0.643, 0.16, 0.02, 0.733, 0.828, 0.901, 0.054, 0.723, 0.816, 0.098, 0.636, 0.345, 0.638, 0.205, 0.033, 0.335, 0.784, 0.015, 0.571, 0.241, 0.9, 0.538, 0.888, 0.49, 0.064, 0.804, 0.541, 0.022, 0.685, 0.339, 0.441, 0.193, 0.086, 0.552, 0.23, 0.674, 0.385, 0.612, 0.31, 0.007, 0.565, 0.758, 0.14, 0.119, 0.946, 0.121, 0.907, 0.358, 0.641, 0.546, 0.175, 0.553, 0.877, 0.802, 0.4, 0.312, 0.377, 0.011, 0.702, 0.324, 0.348, 0.618, 0.962, 0.826, 0.01, 0.737, 0.68, 0.76, 0.027, 0.066, 0.881, 0.552, 0.11, 0.427, 0.322, 0.496, 0.132, 0.907, 0.722, 0.093, 0.754, 0.922, 0.378, 0.277, 0.823, 0.037, 0.266, 0.69, 0.11, 0.273, 0.545, 0.952, 0.088, 0.352, 0.549, 0.747, 0.673, 0.298, 0.02, 0.397, 0.845]
global q = [0.191, 0.999, 0.979, 0.632, 0.999, 0.943, 0.946, 0.74, 0.931, 0.788, 0.399, 0.948, 0.831, 0.752, 0.698, 0.806, 0.684, 0.517, 0.912, 0.757, 0.657, 0.64, 0.622, 0.887, 0.742, 0.952, 0.912, 0.988, 0.829, 0.668, 0.995, 0.461, 0.707, 0.814, 0.876, 0.974, 0.986, 0.876, 0.999, 0.715, 0.366, 0.973, 0.907, 0.985, 0.768, 0.6, 0.697, 0.892, 0.978, 0.968, 0.969, 0.826, 0.839, 0.87, 0.882, 0.876, 0.237, 0.988, 0.531, 0.984, 0.435, 0.154, 0.891, 0.999, 0.72, 0.942, 0.851, 0.737, 0.892, 0.177, 0.875, 0.774, 0.559, 0.807, 0.722, 0.678, 0.604, 0.913, 0.648, 0.725, 0.986, 0.383, 0.974, 0.863, 0.798, 0.766, 0.935, 0.838, 0.982, 0.688, 0.346, 0.823, 0.915, 0.907, 0.848, 0.954, 0.237, 0.97, 0.906, 0.999, 0.865, 0.555, 0.669, 0.757, 0.761, 0.992, 0.991, 0.317, 0.578, 0.356, 0.801, 0.992, 0.41, 0.607, 0.911, 0.982, 0.978, 0.394, 0.863, 0.905, 0.788, 0.134, 0.817, 0.577, 0.963, 0.314, 0.346, 0.837, 0.262, 0.431, 0.933, 0.892, 0.943, 0.436, 0.818, 0.893, 0.734, 0.693, 0.626, 0.972, 0.452, 0.787, 0.593, 0.959, 0.305, 0.882, 0.886, 0.945, 0.777, 0.957, 0.872, 0.1, 0.852, 0.966, 0.308, 0.833, 0.518, 0.639, 0.299, 0.436, 0.728, 0.376, 0.7, 0.752, 0.877, 0.664, 0.069, 0.85, 0.792, 0.579, 0.143, 0.997, 0.624, 0.941, 0.741, 0.948, 0.561, 0.603, 0.877, 0.933, 0.816, 0.452, 0.92, 0.672, 0.63, 0.888, 0.713, 0.599, 0.698, 0.991, 0.907, 0.806, 0.748, 0.99, 0.901, 0.457, 0.967, 0.982, 0.582, 0.712, 0.656, 0.799, 0.579, 0.441, 0.928, 0.836, 0.714, 0.873, 0.956, 0.958, 0.537, 0.932, 0.771, 0.705, 0.909, 0.423, 0.872, 0.623, 0.996, 0.652, 0.584, 0.788, 0.978, 0.712, 0.317, 0.229, 0.714, 0.936]
global origin = 1
global destination = 50