global arcs = [1 10; 1 18; 1 33; 2 8; 2 15; 2 32; 2 47; 3 2; 3 13; 3 29; 3 47; 4 25; 4 44; 4 45; 4 50; 5 19; 5 21; 5 48; 6 4; 6 32; 6 34; 6 37; 7 9; 7 16; 7 20; 7 28; 7 34; 8 10; 8 14; 8 34; 8 38; 8 48; 9 5; 9 35; 9 37; 9 46; 9 50; 10 22; 10 39; 11 12; 11 28; 11 34; 11 37; 11 48; 12 2; 12 3; 12 8; 12 24; 12 34; 12 41; 12 45; 13 8; 13 14; 13 35; 13 36; 13 41; 13 43; 14 24; 14 29; 14 46; 15 17; 15 20; 15 22; 15 29; 15 39; 15 50; 16 11; 16 13; 16 18; 16 34; 16 44; 17 5; 17 11; 17 22; 17 38; 18 7; 19 15; 19 23; 20 8; 20 9; 20 18; 20 25; 20 43; 21 7; 21 13; 21 48; 21 50; 22 8; 22 13; 22 18; 22 23; 22 38; 23 20; 23 32; 23 46; 23 48; 24 6; 24 9; 24 12; 24 19; 24 34; 25 6; 25 7; 25 28; 25 31; 25 48; 26 6; 26 13; 26 22; 27 10; 27 16; 27 19; 27 28; 27 48; 28 17; 28 24; 28 26; 28 37; 29 5; 29 7; 29 12; 29 23; 29 33; 30 6; 30 9; 30 11; 30 31; 30 40; 30 46; 31 5; 31 30; 31 34; 32 2; 32 8; 32 16; 32 25; 32 48; 32 50; 33 7; 33 24; 33 27; 34 13; 34 21; 34 26; 34 38; 34 43; 35 28; 35 31; 35 50; 36 3; 36 5; 36 12; 36 18; 36 21; 36 26; 36 42; 36 46; 36 49; 37 4; 37 11; 37 12; 37 18; 37 44; 37 47; 37 50; 38 29; 38 30; 38 45; 38 46; 39 10; 39 33; 39 38; 39 46; 39 48; 40 6; 40 13; 40 18; 40 27; 40 32; 40 41; 41 15; 41 25; 41 34; 41 38; 42 6; 42 10; 42 16; 42 46; 43 12; 43 22; 43 25; 43 28; 43 40; 43 50; 44 7; 44 12; 44 15; 44 47; 45 25; 45 26; 45 47; 46 9; 46 12; 46 17; 46 25; 46 31; 46 34; 46 36; 47 6; 47 35; 47 36; 48 10; 48 12; 48 32; 48 33; 48 40; 48 45; 49 3; 49 7; 49 18; 49 20; 49 25; 49 32; 49 35; 49 36]
global d_x = [1.0, 2.0, 3.0, 8.0, 2.0, 1.0, 3.0, 4.0, 4.0, 10.0, 1.0, 7.0, 5.0, 6.0, 7.0, 10.0, 3.0, 3.0, 9.0, 3.0, 4.0, 1.0, 3.0, 6.0, 2.0, 10.0, 7.0, 5.0, 7.0, 9.0, 5.0, 9.0, 2.0, 6.0, 4.0, 1.0, 4.0, 1.0, 5.0, 9.0, 10.0, 8.0, 1.0, 8.0, 3.0, 8.0, 6.0, 6.0, 5.0, 2.0, 6.0, 9.0, 9.0, 8.0, 10.0, 4.0, 10.0, 7.0, 4.0, 6.0, 5.0, 6.0, 6.0, 5.0, 1.0, 8.0, 10.0, 1.0, 8.0, 4.0, 7.0, 1.0, 5.0, 1.0, 1.0, 3.0, 3.0, 8.0, 10.0, 10.0, 4.0, 3.0, 4.0, 8.0, 6.0, 9.0, 5.0, 6.0, 8.0, 6.0, 4.0, 1.0, 10.0, 10.0, 8.0, 7.0, 8.0, 3.0, 3.0, 1.0, 1.0, 8.0, 4.0, 7.0, 1.0, 5.0, 5.0, 1.0, 3.0, 10.0, 2.0, 9.0, 3.0, 7.0, 6.0, 10.0, 6.0, 6.0, 9.0, 7.0, 9.0, 6.0, 6.0, 3.0, 1.0, 10.0, 8.0, 6.0, 6.0, 2.0, 1.0, 5.0, 1.0, 6.0, 7.0, 10.0, 3.0, 3.0, 6.0, 5.0, 7.0, 2.0, 5.0, 8.0, 6.0, 5.0, 4.0, 8.0, 3.0, 2.0, 7.0, 9.0, 1.0, 8.0, 8.0, 7.0, 8.0, 1.0, 7.0, 6.0, 5.0, 6.0, 8.0, 5.0, 1.0, 3.0, 1.0, 8.0, 9.0, 2.0, 9.0, 1.0, 5.0, 5.0, 9.0, 7.0, 5.0, 4.0, 4.0, 6.0, 6.0, 1.0, 1.0, 6.0, 5.0, 2.0, 7.0, 7.0, 5.0, 4.0, 9.0, 4.0, 10.0, 9.0, 7.0, 9.0, 8.0, 4.0, 3.0, 5.0, 3.0, 10.0, 4.0, 5.0, 1.0, 9.0, 7.0, 7.0, 2.0, 6.0, 1.0, 3.0, 7.0, 10.0, 5.0, 9.0, 4.0, 3.0, 9.0, 1.0, 2.0, 4.0, 7.0, 5.0, 7.0]
global b_x = 5
global d_y = [8.0, 3.0, 9.0, 8.0, 1.0, 1.0, 3.0, 9.0, 1.0, 3.0, 8.0, 10.0, 2.0, 10.0, 2.0, 8.0, 1.0, 7.0, 3.0, 10.0, 7.0, 4.0, 3.0, 1.0, 3.0, 6.0, 9.0, 1.0, 10.0, 9.0, 3.0, 7.0, 3.0, 7.0, 10.0, 9.0, 1.0, 8.0, 9.0, 4.0, 2.0, 8.0, 10.0, 4.0, 6.0, 10.0, 9.0, 3.0, 5.0, 8.0, 5.0, 3.0, 10.0, 2.0, 7.0, 6.0, 5.0, 1.0, 2.0, 8.0, 1.0, 2.0, 8.0, 2.0, 3.0, 7.0, 9.0, 3.0, 3.0, 9.0, 5.0, 3.0, 4.0, 5.0, 4.0, 8.0, 4.0, 8.0, 4.0, 3.0, 3.0, 10.0, 4.0, 6.0, 4.0, 5.0, 3.0, 2.0, 2.0, 4.0, 2.0, 5.0, 8.0, 6.0, 6.0, 5.0, 8.0, 7.0, 4.0, 3.0, 4.0, 10.0, 6.0, 4.0, 3.0, 1.0, 5.0, 3.0, 8.0, 10.0, 4.0, 9.0, 1.0, 9.0, 1.0, 6.0, 8.0, 9.0, 10.0, 1.0, 7.0, 8.0, 10.0, 9.0, 3.0, 9.0, 7.0, 3.0, 4.0, 5.0, 10.0, 1.0, 8.0, 1.0, 10.0, 5.0, 8.0, 8.0, 7.0, 1.0, 10.0, 3.0, 3.0, 8.0, 8.0, 7.0, 8.0, 2.0, 8.0, 1.0, 4.0, 6.0, 10.0, 5.0, 1.0, 8.0, 6.0, 10.0, 4.0, 1.0, 3.0, 6.0, 1.0, 1.0, 10.0, 9.0, 8.0, 5.0, 10.0, 9.0, 10.0, 7.0, 10.0, 9.0, 8.0, 4.0, 3.0, 7.0, 1.0, 1.0, 1.0, 4.0, 3.0, 3.0, 2.0, 3.0, 4.0, 8.0, 9.0, 4.0, 8.0, 6.0, 1.0, 6.0, 10.0, 6.0, 10.0, 5.0, 3.0, 6.0, 4.0, 3.0, 3.0, 1.0, 9.0, 10.0, 8.0, 2.0, 10.0, 2.0, 8.0, 8.0, 8.0, 1.0, 2.0, 1.0, 6.0, 6.0, 1.0, 2.0, 10.0, 4.0, 1.0, 2.0, 8.0]
global b_y = 10
global p = [0.809, 0.936, 0.449, 0.559, 0.847, 0.368, 0.317, 0.466, 0.607, 0.636, 0.446, 0.997, 0.939, 0.799, 0.835, 0.204, 0.842, 0.149, 0.636, 0.105, 0.321, 0.251, 0.831, 0.498, 0.104, 0.346, 0.108, 0.088, 0.178, 0.295, 0.857, 0.691, 0.963, 0.478, 0.273, 0.733, 0.198, 0.377, 0.856, 0.826, 0.525, 0.67, 0.71, 0.334, 0.246, 0.124, 0.35, 0.822, 0.8, 0.36, 0.81, 0.561, 0.338, 0.857, 0.47, 0.216, 0.729, 0.201, 0.752, 0.006, 0.791, 0.493, 0.363, 0.463, 0.812, 0.489, 0.957, 0.813, 0.57, 0.305, 0.119, 0.847, 0.89, 0.562, 0.458, 0.241, 0.827, 0.722, 0.436, 0.603, 0.403, 0.887, 0.004, 0.269, 0.986, 0.02, 0.704, 0.714, 0.442, 0.372, 0.001, 0.708, 0.105, 0.591, 0.647, 0.921, 0.667, 0.292, 0.147, 0.945, 0.657, 0.62, 0.981, 0.474, 0.733, 0.764, 0.417, 0.728, 0.453, 0.993, 0.047, 0.249, 0.771, 0.942, 0.011, 0.833, 0.991, 0.764, 0.084, 0.111, 0.366, 0.994, 0.461, 0.303, 0.458, 0.885, 0.931, 0.678, 0.067, 0.496, 0.055, 0.793, 0.657, 0.622, 0.832, 0.933, 0.905, 0.411, 0.114, 0.616, 0.332, 0.963, 0.01, 0.911, 0.093, 0.322, 0.542, 0.881, 0.482, 0.039, 0.729, 0.06, 0.775, 0.492, 0.684, 0.096, 0.248, 0.552, 0.001, 0.793, 0.953, 0.223, 0.318, 0.34, 0.187, 0.616, 0.357, 0.057, 0.343, 0.734, 0.047, 0.166, 0.921, 0.841, 0.617, 0.654, 0.167, 0.91, 0.217, 0.376, 0.934, 0.525, 0.515, 0.887, 0.42, 0.319, 0.966, 0.354, 0.988, 0.192, 0.344, 0.794, 0.77, 0.514, 0.176, 0.55, 0.564, 0.896, 0.322, 0.641, 0.115, 0.868, 0.864, 0.88, 0.473, 0.767, 0.607, 0.772, 0.883, 0.208, 0.023, 0.431, 0.506, 0.707, 0.939, 0.75, 0.394, 0.615, 0.606, 0.091, 0.365, 0.779, 0.156, 0.937, 0.579]
global q = [0.903, 0.938, 0.506, 0.842, 0.962, 0.919, 0.896, 0.95, 0.77, 0.946, 0.69, 0.997, 0.948, 0.827, 0.844, 0.716, 0.858, 0.46, 0.81, 0.55, 0.773, 0.916, 0.891, 0.531, 0.155, 0.512, 0.395, 0.361, 0.935, 0.842, 0.922, 0.885, 0.998, 0.556, 0.986, 0.913, 0.77, 0.812, 0.961, 0.866, 0.798, 0.894, 0.791, 0.938, 0.662, 0.459, 0.976, 0.876, 0.968, 0.837, 0.962, 0.619, 0.466, 0.868, 0.678, 0.649, 0.747, 0.899, 0.902, 0.983, 0.795, 0.821, 0.745, 0.785, 0.863, 0.995, 0.961, 0.885, 0.734, 0.816, 0.2, 0.923, 0.938, 0.685, 0.899, 0.833, 0.907, 0.766, 0.687, 0.796, 0.963, 0.923, 0.246, 0.902, 0.988, 0.924, 0.937, 0.889, 0.966, 0.793, 0.945, 0.991, 0.438, 0.736, 0.753, 0.982, 0.834, 0.335, 0.684, 0.959, 0.969, 0.974, 0.991, 0.585, 0.908, 0.929, 0.799, 0.871, 0.837, 0.998, 0.612, 0.631, 0.88, 0.952, 0.351, 0.906, 0.998, 0.979, 0.219, 0.313, 0.55, 0.995, 0.708, 0.776, 0.882, 0.919, 0.967, 0.955, 0.782, 0.496, 0.863, 0.996, 0.907, 0.893, 0.876, 0.991, 0.979, 0.662, 0.3, 0.768, 0.44, 0.993, 0.236, 0.949, 0.3, 0.932, 0.682, 0.984, 0.751, 0.271, 0.798, 0.6, 0.812, 0.774, 0.926, 0.379, 0.758, 0.932, 0.848, 0.952, 0.972, 0.33, 0.602, 0.984, 0.83, 0.765, 0.599, 0.254, 0.43, 0.816, 0.163, 0.915, 0.952, 0.853, 0.664, 0.883, 0.725, 0.941, 0.683, 0.666, 0.981, 0.693, 0.9, 0.981, 0.921, 0.863, 0.992, 0.765, 0.996, 0.738, 0.49, 0.943, 0.87, 0.863, 0.205, 0.596, 0.601, 0.939, 0.453, 0.708, 0.788, 0.919, 0.977, 0.988, 0.961, 0.993, 0.848, 0.937, 0.918, 0.885, 0.299, 0.54, 0.71, 0.911, 0.962, 0.886, 0.527, 0.68, 0.806, 0.745, 0.712, 0.977, 0.399, 0.999, 0.906]
global origin = 1
global destination = 50