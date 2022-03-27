global arcs = [1 5; 1 7; 1 16; 1 25; 1 39; 1 50; 2 14; 2 20; 2 25; 2 35; 2 36; 2 43; 3 24; 3 35; 4 5; 4 13; 4 16; 4 25; 4 36; 5 4; 5 15; 5 19; 5 21; 5 28; 5 30; 5 31; 5 32; 6 8; 6 23; 6 50; 7 10; 7 12; 7 13; 7 19; 7 22; 7 29; 7 35; 7 36; 7 49; 8 3; 8 5; 8 21; 8 25; 8 28; 8 40; 8 45; 9 10; 9 21; 9 37; 9 48; 10 9; 10 21; 11 5; 12 8; 12 9; 12 10; 12 13; 12 22; 12 25; 12 34; 12 50; 13 2; 13 15; 13 21; 13 26; 13 27; 13 29; 13 32; 13 50; 14 3; 14 16; 14 32; 15 21; 15 39; 15 40; 16 7; 16 8; 16 14; 16 19; 16 26; 17 9; 17 15; 17 37; 17 38; 17 45; 18 2; 18 17; 18 30; 19 8; 19 9; 19 10; 19 28; 19 30; 19 38; 20 5; 20 11; 20 37; 21 2; 21 18; 21 34; 21 36; 21 44; 22 12; 22 31; 22 40; 23 17; 23 33; 23 46; 23 49; 24 15; 24 18; 24 21; 24 23; 24 26; 24 30; 24 37; 25 2; 25 4; 25 5; 25 32; 25 36; 26 10; 26 16; 26 27; 26 45; 26 50; 27 10; 27 18; 27 26; 27 30; 27 48; 28 4; 28 10; 28 20; 28 33; 28 46; 29 16; 29 27; 29 41; 30 6; 30 37; 30 40; 31 24; 31 32; 31 33; 31 38; 31 44; 31 49; 32 2; 32 17; 32 19; 32 21; 32 27; 32 30; 32 31; 32 38; 32 39; 33 3; 33 6; 33 13; 33 17; 33 25; 34 5; 34 21; 34 40; 34 43; 35 12; 35 17; 35 24; 35 44; 36 9; 36 13; 36 14; 36 25; 37 4; 37 17; 38 26; 38 34; 38 47; 39 2; 39 18; 39 23; 39 31; 40 9; 40 37; 41 5; 41 10; 41 12; 41 24; 41 35; 41 47; 41 49; 42 31; 42 33; 43 8; 43 12; 43 13; 43 38; 43 45; 44 4; 44 18; 44 26; 44 28; 44 42; 44 43; 44 47; 45 22; 45 25; 45 27; 45 33; 45 34; 45 36; 45 38; 46 15; 46 17; 46 18; 46 24; 46 30; 46 32; 47 14; 47 21; 47 25; 48 10; 48 26; 49 5; 49 31; 49 33; 49 43; 49 46; 49 48]
global d_x = [8.0, 10.0, 1.0, 4.0, 9.0, 6.0, 7.0, 8.0, 2.0, 7.0, 10.0, 9.0, 4.0, 4.0, 7.0, 5.0, 9.0, 4.0, 8.0, 6.0, 3.0, 5.0, 10.0, 7.0, 5.0, 7.0, 4.0, 6.0, 5.0, 5.0, 7.0, 10.0, 6.0, 3.0, 8.0, 4.0, 7.0, 3.0, 5.0, 5.0, 2.0, 7.0, 5.0, 2.0, 5.0, 7.0, 2.0, 6.0, 10.0, 6.0, 6.0, 6.0, 2.0, 9.0, 3.0, 6.0, 4.0, 1.0, 9.0, 6.0, 6.0, 5.0, 5.0, 1.0, 6.0, 7.0, 10.0, 2.0, 6.0, 1.0, 8.0, 1.0, 4.0, 1.0, 7.0, 10.0, 5.0, 2.0, 8.0, 5.0, 6.0, 2.0, 5.0, 2.0, 5.0, 8.0, 4.0, 7.0, 10.0, 5.0, 1.0, 6.0, 7.0, 3.0, 7.0, 4.0, 5.0, 1.0, 6.0, 4.0, 8.0, 2.0, 1.0, 6.0, 4.0, 6.0, 3.0, 4.0, 3.0, 6.0, 9.0, 3.0, 3.0, 5.0, 2.0, 2.0, 8.0, 4.0, 1.0, 9.0, 3.0, 5.0, 8.0, 7.0, 6.0, 10.0, 5.0, 10.0, 3.0, 10.0, 1.0, 5.0, 6.0, 10.0, 9.0, 6.0, 7.0, 6.0, 2.0, 1.0, 8.0, 9.0, 2.0, 10.0, 2.0, 5.0, 7.0, 5.0, 6.0, 2.0, 10.0, 9.0, 9.0, 9.0, 2.0, 1.0, 8.0, 5.0, 7.0, 3.0, 5.0, 6.0, 9.0, 7.0, 10.0, 2.0, 6.0, 2.0, 2.0, 2.0, 4.0, 3.0, 3.0, 7.0, 7.0, 8.0, 2.0, 1.0, 4.0, 1.0, 4.0, 5.0, 7.0, 6.0, 8.0, 10.0, 6.0, 3.0, 6.0, 5.0, 4.0, 6.0, 3.0, 8.0, 8.0, 10.0, 2.0, 8.0, 5.0, 10.0, 8.0, 5.0, 10.0, 10.0, 3.0, 4.0, 2.0, 8.0, 6.0, 6.0, 5.0, 4.0, 8.0, 1.0, 2.0, 9.0, 5.0, 2.0, 5.0, 1.0, 7.0, 6.0, 7.0, 10.0, 7.0, 8.0, 5.0, 3.0, 4.0, 1.0]
global b_x = 5
global d_y = [10.0, 8.0, 9.0, 5.0, 6.0, 2.0, 4.0, 2.0, 2.0, 7.0, 9.0, 2.0, 1.0, 9.0, 3.0, 3.0, 10.0, 6.0, 10.0, 10.0, 8.0, 8.0, 6.0, 10.0, 7.0, 3.0, 7.0, 6.0, 7.0, 2.0, 8.0, 4.0, 4.0, 2.0, 7.0, 6.0, 10.0, 3.0, 9.0, 8.0, 10.0, 3.0, 7.0, 6.0, 3.0, 4.0, 9.0, 1.0, 1.0, 2.0, 2.0, 8.0, 9.0, 6.0, 7.0, 10.0, 3.0, 1.0, 8.0, 1.0, 4.0, 4.0, 8.0, 9.0, 5.0, 8.0, 8.0, 7.0, 2.0, 8.0, 6.0, 7.0, 10.0, 10.0, 7.0, 5.0, 4.0, 6.0, 1.0, 2.0, 4.0, 4.0, 6.0, 1.0, 7.0, 3.0, 5.0, 10.0, 4.0, 8.0, 8.0, 5.0, 1.0, 10.0, 5.0, 5.0, 3.0, 2.0, 10.0, 1.0, 3.0, 7.0, 7.0, 1.0, 7.0, 8.0, 4.0, 10.0, 7.0, 4.0, 3.0, 6.0, 2.0, 3.0, 7.0, 2.0, 2.0, 7.0, 4.0, 1.0, 6.0, 3.0, 6.0, 6.0, 8.0, 3.0, 4.0, 8.0, 1.0, 9.0, 7.0, 4.0, 1.0, 6.0, 8.0, 9.0, 3.0, 4.0, 5.0, 10.0, 7.0, 8.0, 8.0, 2.0, 1.0, 4.0, 9.0, 9.0, 6.0, 7.0, 9.0, 1.0, 10.0, 6.0, 1.0, 6.0, 7.0, 6.0, 1.0, 4.0, 2.0, 4.0, 8.0, 3.0, 1.0, 6.0, 8.0, 8.0, 8.0, 4.0, 1.0, 9.0, 8.0, 1.0, 10.0, 7.0, 10.0, 1.0, 3.0, 8.0, 2.0, 1.0, 6.0, 4.0, 6.0, 5.0, 2.0, 9.0, 5.0, 8.0, 8.0, 2.0, 5.0, 10.0, 2.0, 3.0, 7.0, 10.0, 8.0, 4.0, 2.0, 6.0, 6.0, 1.0, 7.0, 2.0, 6.0, 5.0, 7.0, 2.0, 1.0, 10.0, 5.0, 9.0, 5.0, 1.0, 2.0, 1.0, 8.0, 7.0, 3.0, 10.0, 4.0, 10.0, 6.0, 1.0, 2.0, 2.0, 10.0, 5.0]
global b_y = 10
global p = [0.46, 0.912, 0.954, 0.168, 0.28, 0.704, 0.249, 0.819, 0.511, 0.954, 0.512, 0.683, 0.744, 0.838, 0.365, 0.658, 0.474, 0.941, 0.165, 0.316, 0.618, 0.93, 0.599, 0.678, 0.912, 0.961, 0.237, 0.72, 0.839, 0.329, 0.943, 0.953, 0.542, 0.057, 0.056, 0.59, 0.312, 0.67, 0.4, 0.424, 0.673, 0.979, 0.019, 0.848, 0.552, 0.6, 0.103, 0.856, 0.483, 0.788, 0.704, 0.198, 0.637, 0.799, 0.373, 0.205, 0.042, 0.727, 0.678, 0.732, 0.047, 0.495, 0.099, 0.069, 0.633, 0.011, 0.223, 0.41, 0.264, 0.462, 0.953, 0.989, 0.159, 0.766, 0.172, 0.887, 0.506, 0.232, 0.845, 0.95, 0.097, 0.153, 0.016, 0.397, 0.603, 0.955, 0.901, 0.276, 0.316, 0.503, 0.937, 0.901, 0.29, 0.901, 0.445, 0.428, 0.079, 0.64, 0.159, 0.965, 0.04, 0.075, 0.696, 0.153, 0.944, 0.348, 0.799, 0.54, 0.058, 0.352, 0.197, 0.781, 0.388, 0.319, 0.154, 0.298, 0.161, 0.558, 0.443, 0.16, 0.063, 0.475, 0.762, 0.497, 0.875, 0.53, 0.132, 0.817, 0.147, 0.773, 0.096, 0.352, 0.638, 0.553, 0.372, 0.895, 0.341, 0.981, 0.141, 0.461, 0.017, 0.396, 0.832, 0.47, 0.825, 0.1, 0.562, 0.694, 0.908, 0.519, 0.188, 0.409, 0.87, 0.462, 0.164, 0.41, 0.818, 0.678, 0.768, 0.069, 0.534, 0.202, 0.67, 0.749, 0.694, 0.61, 0.279, 0.287, 0.106, 0.334, 0.606, 0.475, 0.834, 0.692, 0.546, 0.873, 0.301, 0.949, 0.378, 0.564, 0.321, 0.969, 0.125, 0.918, 0.213, 0.854, 0.766, 0.977, 0.71, 0.021, 0.6, 0.629, 0.651, 0.45, 0.648, 0.027, 0.643, 0.866, 0.712, 0.08, 0.77, 0.719, 0.794, 0.245, 0.751, 0.011, 0.619, 0.827, 0.203, 0.521, 0.996, 0.743, 0.096, 0.149, 0.711, 0.108, 0.43, 0.164, 0.709, 0.422, 0.868, 0.751, 0.974, 0.404, 0.082, 0.922, 0.343, 0.955, 0.361, 0.053]
global q = [0.508, 0.914, 0.966, 0.884, 0.67, 0.751, 0.6, 0.89, 0.683, 0.997, 0.563, 0.798, 0.787, 0.966, 0.882, 0.859, 0.761, 0.949, 0.989, 0.902, 0.889, 0.968, 0.771, 0.693, 0.937, 0.99, 0.302, 0.953, 0.955, 0.641, 0.975, 0.959, 0.897, 0.967, 0.381, 0.893, 0.756, 0.808, 0.467, 0.49, 0.986, 0.989, 0.556, 0.899, 0.634, 0.967, 0.535, 0.914, 0.82, 0.839, 0.733, 0.273, 0.785, 0.974, 0.71, 0.99, 0.531, 0.894, 0.968, 0.832, 0.631, 0.845, 0.284, 0.549, 0.635, 0.788, 0.617, 0.682, 0.911, 0.788, 0.966, 0.992, 0.69, 0.854, 0.557, 0.96, 0.674, 0.368, 0.911, 0.981, 0.82, 0.495, 0.193, 0.951, 0.974, 0.989, 0.995, 0.693, 0.702, 0.888, 0.973, 0.953, 0.367, 0.946, 0.634, 0.824, 0.463, 0.677, 0.989, 0.996, 0.228, 0.644, 0.901, 0.247, 0.981, 0.508, 0.923, 0.685, 0.478, 0.857, 0.569, 0.883, 0.544, 0.429, 0.183, 0.429, 0.262, 0.905, 0.464, 0.96, 0.988, 0.767, 0.927, 0.649, 0.95, 0.695, 0.413, 0.961, 0.544, 0.867, 0.168, 0.591, 0.836, 0.767, 0.396, 0.919, 0.741, 0.998, 0.782, 0.567, 0.913, 0.651, 0.92, 0.508, 0.845, 0.197, 0.652, 0.765, 0.947, 0.713, 0.231, 0.477, 0.989, 0.827, 0.837, 0.693, 0.936, 0.824, 0.769, 0.846, 0.918, 0.423, 0.991, 0.85, 0.851, 0.837, 0.326, 0.968, 0.208, 0.872, 0.972, 0.814, 0.853, 0.89, 0.621, 0.908, 0.987, 0.961, 0.889, 0.911, 0.548, 0.992, 0.638, 0.993, 0.952, 0.881, 0.809, 0.981, 0.863, 0.686, 0.723, 0.985, 0.83, 0.587, 0.852, 0.198, 0.912, 0.883, 0.821, 0.141, 0.913, 0.888, 0.801, 0.454, 0.818, 0.712, 0.673, 0.85, 0.841, 0.558, 0.999, 0.771, 0.11, 0.833, 0.824, 0.9, 0.555, 0.765, 0.818, 0.768, 0.879, 0.823, 0.982, 0.729, 0.564, 0.991, 0.923, 0.977, 0.362, 0.458]
global origin = 1
global destination = 50