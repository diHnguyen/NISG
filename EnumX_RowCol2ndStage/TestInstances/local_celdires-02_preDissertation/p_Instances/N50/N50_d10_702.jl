global arcs = [1 9; 1 43; 1 45; 2 3; 2 12; 2 23; 2 45; 2 48; 3 10; 3 11; 3 19; 3 20; 3 32; 3 45; 3 49; 4 24; 4 44; 5 2; 5 29; 5 37; 5 45; 6 10; 6 18; 6 37; 7 25; 8 5; 8 10; 8 32; 8 35; 8 37; 8 45; 8 50; 9 31; 9 39; 10 9; 10 14; 10 16; 10 19; 10 39; 11 4; 11 10; 11 14; 11 22; 11 24; 11 36; 12 2; 12 17; 12 29; 12 39; 12 40; 12 43; 13 5; 13 12; 13 33; 13 34; 13 47; 13 49; 14 7; 14 17; 14 47; 14 50; 15 4; 15 6; 15 24; 15 43; 15 49; 16 4; 16 7; 16 11; 16 25; 16 49; 17 7; 17 9; 17 20; 17 43; 18 11; 18 37; 18 39; 18 49; 19 4; 19 13; 19 26; 19 41; 20 11; 20 25; 20 39; 20 44; 21 11; 21 28; 21 31; 21 39; 22 10; 22 36; 22 38; 23 6; 23 10; 23 13; 23 19; 23 27; 23 33; 23 36; 23 38; 24 5; 24 16; 24 21; 24 37; 24 42; 25 15; 25 20; 25 21; 25 23; 26 7; 26 18; 26 23; 26 30; 26 32; 26 37; 27 6; 27 14; 27 19; 27 28; 27 35; 28 16; 29 4; 29 5; 29 11; 29 16; 29 21; 29 23; 29 33; 29 34; 30 32; 30 44; 31 3; 31 6; 31 10; 31 12; 31 19; 31 22; 31 27; 31 34; 31 44; 32 14; 32 41; 32 47; 33 8; 33 9; 33 13; 33 22; 33 23; 34 16; 34 21; 34 22; 34 30; 34 40; 34 41; 34 44; 35 4; 35 11; 35 15; 35 22; 35 23; 35 28; 35 39; 35 41; 36 9; 36 10; 36 21; 36 41; 37 3; 37 5; 37 9; 37 25; 37 27; 38 2; 38 32; 38 40; 38 47; 39 3; 39 13; 39 24; 39 28; 39 36; 40 31; 40 32; 40 44; 41 23; 41 27; 41 30; 42 5; 42 13; 42 28; 42 37; 42 46; 43 19; 43 35; 43 36; 43 44; 44 20; 44 34; 44 37; 45 29; 45 40; 46 6; 46 8; 46 16; 46 28; 46 38; 46 40; 46 50; 47 2; 47 3; 47 5; 47 28; 47 38; 47 42; 48 5; 48 12; 48 14; 48 17; 48 24; 48 30; 48 31; 48 46; 48 49; 48 50; 49 5; 49 24; 49 32; 49 50]
global d_x = [10.0, 4.0, 1.0, 3.0, 2.0, 10.0, 8.0, 7.0, 2.0, 1.0, 1.0, 9.0, 10.0, 9.0, 7.0, 5.0, 2.0, 5.0, 4.0, 9.0, 8.0, 2.0, 1.0, 9.0, 8.0, 1.0, 8.0, 7.0, 3.0, 4.0, 9.0, 1.0, 8.0, 6.0, 5.0, 10.0, 8.0, 2.0, 6.0, 7.0, 3.0, 10.0, 2.0, 4.0, 1.0, 7.0, 3.0, 9.0, 10.0, 7.0, 2.0, 2.0, 4.0, 7.0, 10.0, 10.0, 6.0, 5.0, 5.0, 9.0, 1.0, 8.0, 8.0, 9.0, 6.0, 5.0, 5.0, 9.0, 1.0, 9.0, 2.0, 8.0, 7.0, 6.0, 1.0, 2.0, 1.0, 6.0, 5.0, 10.0, 5.0, 9.0, 7.0, 5.0, 5.0, 9.0, 4.0, 4.0, 7.0, 2.0, 4.0, 6.0, 4.0, 2.0, 6.0, 6.0, 10.0, 4.0, 5.0, 5.0, 3.0, 8.0, 8.0, 10.0, 4.0, 3.0, 2.0, 2.0, 3.0, 10.0, 5.0, 3.0, 5.0, 8.0, 8.0, 6.0, 9.0, 3.0, 5.0, 4.0, 7.0, 4.0, 6.0, 2.0, 4.0, 2.0, 10.0, 10.0, 4.0, 3.0, 6.0, 6.0, 3.0, 4.0, 1.0, 8.0, 6.0, 5.0, 10.0, 5.0, 3.0, 5.0, 9.0, 6.0, 6.0, 5.0, 7.0, 8.0, 10.0, 8.0, 7.0, 7.0, 4.0, 10.0, 7.0, 1.0, 8.0, 7.0, 9.0, 4.0, 5.0, 5.0, 4.0, 3.0, 10.0, 6.0, 10.0, 8.0, 8.0, 3.0, 3.0, 10.0, 9.0, 7.0, 5.0, 4.0, 6.0, 7.0, 7.0, 9.0, 5.0, 9.0, 4.0, 6.0, 7.0, 2.0, 8.0, 3.0, 5.0, 10.0, 9.0, 5.0, 6.0, 1.0, 5.0, 2.0, 7.0, 3.0, 1.0, 8.0, 10.0, 7.0, 8.0, 4.0, 9.0, 9.0, 10.0, 7.0, 10.0, 1.0, 2.0, 8.0, 2.0, 9.0, 8.0, 2.0, 7.0, 2.0, 10.0, 7.0, 7.0, 5.0, 7.0, 10.0, 7.0, 1.0, 3.0, 6.0, 3.0, 9.0]
global b_x = 5
global d_y = [2.0, 10.0, 10.0, 4.0, 5.0, 2.0, 10.0, 3.0, 3.0, 4.0, 6.0, 5.0, 4.0, 6.0, 5.0, 8.0, 3.0, 1.0, 7.0, 8.0, 8.0, 2.0, 3.0, 9.0, 7.0, 3.0, 7.0, 8.0, 8.0, 7.0, 1.0, 5.0, 9.0, 8.0, 2.0, 1.0, 5.0, 7.0, 8.0, 3.0, 4.0, 10.0, 3.0, 2.0, 9.0, 10.0, 1.0, 7.0, 6.0, 1.0, 10.0, 4.0, 5.0, 6.0, 2.0, 2.0, 7.0, 6.0, 2.0, 6.0, 5.0, 9.0, 9.0, 10.0, 1.0, 8.0, 10.0, 9.0, 8.0, 6.0, 2.0, 1.0, 2.0, 6.0, 9.0, 10.0, 2.0, 8.0, 4.0, 9.0, 6.0, 10.0, 2.0, 5.0, 7.0, 7.0, 1.0, 9.0, 2.0, 2.0, 6.0, 1.0, 2.0, 2.0, 1.0, 7.0, 9.0, 2.0, 6.0, 10.0, 1.0, 1.0, 4.0, 8.0, 2.0, 1.0, 5.0, 9.0, 1.0, 7.0, 9.0, 4.0, 7.0, 9.0, 2.0, 2.0, 5.0, 7.0, 10.0, 8.0, 2.0, 4.0, 7.0, 3.0, 3.0, 9.0, 3.0, 7.0, 5.0, 1.0, 8.0, 8.0, 9.0, 2.0, 1.0, 8.0, 4.0, 2.0, 6.0, 3.0, 9.0, 1.0, 4.0, 8.0, 8.0, 4.0, 3.0, 4.0, 9.0, 1.0, 10.0, 8.0, 3.0, 1.0, 3.0, 2.0, 2.0, 9.0, 2.0, 5.0, 4.0, 5.0, 4.0, 3.0, 1.0, 7.0, 9.0, 1.0, 4.0, 6.0, 7.0, 10.0, 9.0, 5.0, 9.0, 3.0, 7.0, 7.0, 10.0, 4.0, 2.0, 1.0, 9.0, 5.0, 3.0, 8.0, 7.0, 5.0, 9.0, 3.0, 2.0, 7.0, 2.0, 2.0, 3.0, 3.0, 3.0, 5.0, 3.0, 4.0, 6.0, 1.0, 1.0, 1.0, 7.0, 5.0, 7.0, 2.0, 3.0, 1.0, 10.0, 2.0, 5.0, 4.0, 2.0, 2.0, 10.0, 2.0, 3.0, 10.0, 5.0, 8.0, 1.0, 9.0, 4.0, 6.0, 8.0, 9.0, 2.0, 8.0]
global b_y = 10
global p = [0.658, 0.397, 0.001, 0.228, 0.395, 0.033, 0.601, 0.383, 0.931, 0.875, 0.939, 0.225, 0.153, 0.335, 0.55, 0.92, 0.756, 0.129, 0.111, 0.929, 0.58, 0.744, 0.397, 0.376, 0.758, 0.212, 0.26, 0.199, 0.737, 0.82, 0.016, 0.881, 0.023, 0.788, 0.25, 0.21, 0.363, 0.725, 0.258, 0.848, 0.419, 0.639, 0.656, 0.32, 0.266, 0.919, 0.633, 0.342, 0.233, 0.903, 0.987, 0.04, 0.218, 0.755, 0.968, 0.649, 0.186, 0.639, 0.868, 0.152, 0.023, 0.323, 0.257, 0.967, 0.563, 0.196, 0.401, 0.929, 0.078, 0.981, 0.659, 0.829, 0.442, 0.225, 0.523, 0.377, 0.286, 0.045, 0.269, 0.045, 0.563, 0.502, 0.876, 0.457, 0.567, 0.616, 0.921, 0.836, 0.864, 0.496, 0.001, 0.489, 0.204, 0.387, 0.381, 0.267, 0.645, 0.907, 0.436, 0.878, 0.13, 0.33, 0.886, 0.166, 0.521, 0.673, 0.42, 0.602, 0.54, 0.327, 0.616, 0.632, 0.702, 0.961, 0.006, 0.945, 0.589, 0.292, 0.18, 0.421, 0.442, 0.421, 0.054, 0.069, 0.405, 0.087, 0.389, 0.421, 0.186, 0.552, 0.755, 0.743, 0.213, 0.51, 0.3, 0.139, 0.86, 0.127, 0.69, 0.438, 0.565, 0.137, 0.909, 0.604, 0.354, 0.973, 0.258, 0.698, 0.262, 0.726, 0.682, 0.067, 0.044, 0.592, 0.706, 0.555, 0.601, 0.308, 0.751, 0.49, 0.881, 0.826, 0.768, 0.919, 0.208, 0.062, 0.245, 0.658, 0.178, 0.698, 0.921, 0.182, 0.078, 0.65, 0.313, 0.207, 0.229, 0.48, 0.376, 0.524, 0.479, 0.865, 0.816, 0.155, 0.957, 0.089, 0.169, 0.028, 0.157, 0.37, 0.374, 0.48, 0.699, 0.199, 0.682, 0.924, 0.458, 0.144, 0.975, 0.52, 0.186, 0.878, 0.457, 0.077, 0.433, 0.738, 0.259, 0.432, 0.054, 0.781, 0.232, 0.965, 0.956, 0.461, 0.839, 0.965, 0.708, 0.163, 0.847, 0.366, 0.666, 0.336, 0.987, 0.273, 0.878, 0.006, 0.149, 0.899, 0.865, 0.899]
global q = [0.873, 0.777, 0.516, 0.268, 0.526, 0.811, 0.871, 0.492, 0.933, 0.981, 0.94, 0.997, 0.813, 0.836, 0.926, 0.929, 0.828, 0.717, 0.791, 0.938, 0.981, 0.773, 0.953, 0.712, 0.983, 0.771, 0.689, 0.837, 0.817, 0.842, 0.951, 0.978, 0.26, 0.948, 0.62, 0.979, 0.381, 0.783, 0.93, 0.93, 0.946, 0.912, 0.704, 0.76, 0.368, 0.94, 0.803, 0.619, 0.328, 0.958, 0.987, 0.689, 0.7, 0.976, 0.974, 0.713, 0.73, 0.745, 0.934, 0.871, 0.555, 0.88, 0.679, 0.996, 0.99, 0.476, 0.635, 0.942, 0.852, 0.984, 0.991, 0.935, 0.51, 0.617, 0.676, 0.394, 0.443, 0.299, 0.777, 0.934, 0.791, 0.833, 0.922, 0.975, 0.63, 0.716, 0.966, 0.851, 0.916, 0.603, 0.878, 0.582, 0.801, 0.401, 0.769, 0.387, 0.692, 0.939, 0.86, 0.917, 0.316, 0.679, 0.93, 0.402, 0.851, 0.892, 0.572, 0.692, 0.688, 0.629, 0.769, 0.908, 0.822, 0.994, 0.932, 0.991, 0.604, 0.735, 0.947, 0.578, 0.744, 0.73, 0.152, 0.502, 0.936, 0.113, 0.872, 0.814, 0.457, 0.97, 0.794, 0.943, 0.783, 0.651, 0.559, 0.703, 0.869, 0.648, 0.975, 0.988, 0.772, 0.95, 0.926, 0.887, 0.536, 0.977, 0.279, 0.881, 0.548, 0.749, 0.777, 0.542, 0.37, 0.671, 0.838, 0.749, 0.938, 0.584, 0.808, 0.877, 0.995, 0.849, 0.776, 0.937, 0.928, 0.241, 0.476, 0.983, 0.445, 0.98, 0.938, 0.652, 0.483, 0.888, 0.779, 0.89, 0.293, 0.652, 0.475, 0.734, 0.598, 0.873, 0.889, 0.549, 0.995, 0.153, 0.792, 0.813, 0.911, 0.939, 0.863, 0.688, 0.802, 0.87, 0.958, 0.934, 0.895, 0.627, 0.993, 0.899, 0.625, 0.979, 0.953, 0.896, 0.646, 0.868, 0.317, 0.569, 0.963, 0.838, 0.996, 0.984, 0.964, 0.482, 0.882, 0.965, 0.939, 0.427, 0.919, 0.769, 0.809, 0.753, 0.989, 0.441, 0.918, 0.679, 0.195, 0.901, 0.943, 0.969]
global origin = 1
global destination = 50