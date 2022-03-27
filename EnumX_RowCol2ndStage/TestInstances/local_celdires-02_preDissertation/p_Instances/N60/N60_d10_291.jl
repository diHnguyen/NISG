global arcs = [1 20; 1 21; 1 23; 1 25; 1 37; 1 56; 1 59; 2 17; 2 37; 2 39; 2 42; 2 46; 2 48; 2 49; 2 51; 2 57; 3 5; 3 7; 3 29; 4 40; 4 46; 4 59; 5 3; 5 11; 5 22; 5 33; 5 35; 5 53; 5 56; 6 10; 6 11; 6 17; 6 35; 6 43; 7 54; 8 4; 8 16; 8 55; 8 58; 8 60; 9 5; 9 19; 9 22; 9 29; 9 34; 9 46; 9 54; 10 3; 10 15; 10 23; 10 44; 10 49; 10 51; 11 5; 11 7; 11 16; 11 18; 11 28; 11 29; 11 42; 11 48; 11 52; 11 59; 12 3; 12 47; 12 56; 13 11; 13 18; 13 36; 13 44; 13 50; 14 9; 14 11; 14 17; 14 19; 14 20; 14 22; 14 24; 14 27; 14 44; 15 12; 15 16; 15 44; 16 3; 16 12; 16 15; 16 18; 16 32; 16 37; 16 43; 16 49; 17 5; 17 14; 17 25; 17 32; 17 35; 17 41; 17 43; 17 56; 18 13; 18 24; 18 29; 18 55; 18 58; 19 29; 19 33; 19 47; 19 58; 20 5; 20 8; 20 16; 20 17; 20 22; 20 34; 20 40; 20 46; 20 47; 20 51; 20 56; 21 14; 21 30; 21 42; 21 52; 22 5; 22 26; 22 47; 22 49; 22 54; 23 18; 23 25; 23 53; 23 55; 23 58; 24 8; 24 11; 24 29; 24 46; 24 48; 24 54; 25 7; 25 16; 25 22; 25 31; 25 38; 26 22; 26 36; 27 2; 27 12; 27 15; 27 38; 27 46; 28 7; 28 8; 28 31; 28 59; 29 7; 29 10; 29 23; 29 35; 29 38; 29 60; 30 2; 30 18; 30 46; 30 50; 30 56; 31 20; 31 30; 31 42; 31 46; 32 5; 32 15; 32 27; 32 37; 32 48; 32 54; 32 57; 33 15; 33 30; 33 38; 33 48; 33 59; 34 8; 34 9; 34 12; 34 20; 34 24; 34 26; 34 27; 34 41; 34 46; 35 25; 35 30; 35 45; 35 46; 35 53; 36 6; 36 9; 36 24; 36 27; 36 40; 36 47; 36 58; 37 39; 37 44; 37 48; 37 57; 38 11; 38 27; 38 30; 38 56; 38 57; 39 4; 39 9; 39 19; 39 35; 39 48; 40 4; 40 11; 40 19; 40 27; 40 28; 41 11; 41 16; 41 24; 41 26; 41 31; 41 49; 41 58; 42 5; 42 34; 42 57; 43 9; 43 35; 43 46; 43 53; 44 8; 44 10; 45 7; 45 13; 45 22; 45 26; 45 28; 46 7; 46 8; 46 17; 46 36; 46 50; 46 51; 47 16; 47 20; 47 26; 47 29; 47 36; 47 39; 47 50; 47 54; 48 52; 49 5; 49 14; 49 35; 49 36; 49 38; 49 52; 50 10; 50 12; 50 14; 50 23; 50 39; 50 55; 50 57; 51 3; 51 52; 51 58; 52 5; 52 7; 52 42; 53 9; 53 39; 54 18; 54 21; 54 22; 54 25; 54 34; 54 42; 54 43; 54 47; 54 56; 55 3; 55 14; 55 20; 55 32; 55 58; 55 59; 56 6; 56 20; 56 22; 56 24; 56 28; 56 31; 56 34; 56 45; 56 46; 57 8; 57 13; 57 17; 57 51; 57 52; 58 3; 58 10; 58 17; 58 28; 58 44; 58 51; 59 2; 59 18; 59 48; 59 60]
global d_x = [1.0, 9.0, 4.0, 10.0, 9.0, 4.0, 3.0, 3.0, 1.0, 6.0, 5.0, 2.0, 10.0, 4.0, 4.0, 2.0, 1.0, 2.0, 2.0, 10.0, 8.0, 4.0, 2.0, 8.0, 6.0, 1.0, 2.0, 10.0, 4.0, 5.0, 1.0, 4.0, 9.0, 7.0, 4.0, 10.0, 10.0, 4.0, 3.0, 1.0, 9.0, 5.0, 3.0, 4.0, 3.0, 3.0, 2.0, 6.0, 1.0, 8.0, 2.0, 7.0, 6.0, 9.0, 9.0, 9.0, 4.0, 3.0, 1.0, 9.0, 6.0, 6.0, 10.0, 2.0, 5.0, 6.0, 5.0, 2.0, 2.0, 3.0, 3.0, 6.0, 8.0, 5.0, 7.0, 5.0, 4.0, 5.0, 8.0, 5.0, 9.0, 9.0, 7.0, 9.0, 1.0, 8.0, 10.0, 8.0, 1.0, 7.0, 6.0, 9.0, 5.0, 6.0, 2.0, 4.0, 8.0, 5.0, 2.0, 1.0, 7.0, 4.0, 4.0, 5.0, 7.0, 3.0, 6.0, 10.0, 7.0, 10.0, 7.0, 5.0, 1.0, 10.0, 10.0, 1.0, 10.0, 1.0, 9.0, 2.0, 3.0, 5.0, 8.0, 3.0, 2.0, 10.0, 3.0, 5.0, 1.0, 3.0, 6.0, 1.0, 5.0, 3.0, 8.0, 5.0, 7.0, 3.0, 7.0, 1.0, 8.0, 6.0, 5.0, 7.0, 2.0, 6.0, 7.0, 5.0, 10.0, 4.0, 4.0, 2.0, 6.0, 6.0, 6.0, 5.0, 1.0, 1.0, 6.0, 2.0, 1.0, 7.0, 9.0, 10.0, 5.0, 9.0, 4.0, 2.0, 10.0, 6.0, 4.0, 5.0, 3.0, 8.0, 1.0, 4.0, 10.0, 3.0, 10.0, 8.0, 3.0, 3.0, 9.0, 7.0, 9.0, 2.0, 6.0, 1.0, 7.0, 3.0, 3.0, 1.0, 6.0, 9.0, 9.0, 4.0, 1.0, 2.0, 8.0, 2.0, 6.0, 9.0, 10.0, 8.0, 10.0, 5.0, 8.0, 2.0, 7.0, 5.0, 1.0, 7.0, 1.0, 6.0, 7.0, 7.0, 7.0, 5.0, 8.0, 10.0, 5.0, 10.0, 3.0, 5.0, 4.0, 3.0, 6.0, 5.0, 4.0, 4.0, 3.0, 8.0, 7.0, 8.0, 7.0, 4.0, 8.0, 10.0, 3.0, 9.0, 8.0, 9.0, 3.0, 2.0, 9.0, 4.0, 9.0, 8.0, 4.0, 7.0, 4.0, 10.0, 3.0, 7.0, 4.0, 6.0, 5.0, 9.0, 1.0, 7.0, 3.0, 3.0, 5.0, 2.0, 9.0, 7.0, 6.0, 3.0, 9.0, 1.0, 8.0, 3.0, 9.0, 9.0, 4.0, 6.0, 9.0, 10.0, 3.0, 8.0, 7.0, 8.0, 7.0, 1.0, 2.0, 8.0, 2.0, 8.0, 7.0, 1.0, 3.0, 3.0, 10.0, 3.0, 8.0, 5.0, 7.0, 1.0, 6.0, 2.0, 6.0, 2.0, 3.0, 8.0, 8.0, 9.0, 1.0, 8.0, 2.0, 6.0, 6.0, 10.0, 9.0, 3.0, 1.0, 8.0, 1.0, 6.0]
global b_x = 5
global d_y = [4.0, 2.0, 10.0, 6.0, 3.0, 2.0, 3.0, 5.0, 4.0, 10.0, 4.0, 5.0, 4.0, 7.0, 9.0, 4.0, 4.0, 9.0, 6.0, 7.0, 4.0, 10.0, 2.0, 2.0, 4.0, 9.0, 3.0, 6.0, 10.0, 4.0, 10.0, 2.0, 3.0, 9.0, 4.0, 3.0, 5.0, 10.0, 6.0, 10.0, 8.0, 4.0, 9.0, 5.0, 6.0, 6.0, 8.0, 9.0, 4.0, 3.0, 7.0, 9.0, 7.0, 2.0, 7.0, 6.0, 3.0, 7.0, 10.0, 5.0, 2.0, 7.0, 9.0, 1.0, 10.0, 1.0, 7.0, 8.0, 6.0, 4.0, 2.0, 4.0, 6.0, 3.0, 2.0, 4.0, 10.0, 1.0, 4.0, 6.0, 2.0, 10.0, 9.0, 8.0, 6.0, 6.0, 9.0, 10.0, 3.0, 9.0, 6.0, 10.0, 9.0, 1.0, 6.0, 10.0, 9.0, 3.0, 5.0, 4.0, 9.0, 1.0, 4.0, 1.0, 7.0, 5.0, 1.0, 4.0, 6.0, 10.0, 3.0, 8.0, 8.0, 1.0, 3.0, 5.0, 4.0, 4.0, 5.0, 4.0, 5.0, 2.0, 9.0, 5.0, 9.0, 5.0, 2.0, 10.0, 4.0, 1.0, 7.0, 6.0, 4.0, 2.0, 9.0, 9.0, 9.0, 6.0, 4.0, 8.0, 5.0, 4.0, 4.0, 9.0, 1.0, 7.0, 3.0, 5.0, 10.0, 2.0, 8.0, 9.0, 4.0, 5.0, 6.0, 6.0, 6.0, 3.0, 9.0, 2.0, 9.0, 6.0, 2.0, 5.0, 5.0, 7.0, 3.0, 8.0, 6.0, 9.0, 5.0, 7.0, 1.0, 4.0, 8.0, 9.0, 5.0, 1.0, 3.0, 4.0, 6.0, 8.0, 6.0, 9.0, 10.0, 4.0, 9.0, 5.0, 1.0, 1.0, 8.0, 3.0, 4.0, 9.0, 9.0, 9.0, 1.0, 6.0, 7.0, 2.0, 3.0, 4.0, 7.0, 6.0, 5.0, 1.0, 7.0, 3.0, 8.0, 4.0, 10.0, 9.0, 6.0, 3.0, 10.0, 1.0, 9.0, 3.0, 8.0, 7.0, 4.0, 9.0, 4.0, 10.0, 1.0, 6.0, 5.0, 9.0, 7.0, 10.0, 2.0, 2.0, 3.0, 3.0, 6.0, 2.0, 2.0, 5.0, 7.0, 1.0, 3.0, 10.0, 3.0, 6.0, 6.0, 7.0, 2.0, 10.0, 10.0, 3.0, 10.0, 4.0, 6.0, 2.0, 2.0, 2.0, 9.0, 2.0, 8.0, 2.0, 4.0, 10.0, 8.0, 5.0, 5.0, 8.0, 7.0, 5.0, 6.0, 3.0, 7.0, 4.0, 4.0, 1.0, 4.0, 3.0, 7.0, 4.0, 4.0, 7.0, 6.0, 8.0, 6.0, 6.0, 10.0, 7.0, 4.0, 3.0, 5.0, 5.0, 9.0, 7.0, 7.0, 3.0, 5.0, 8.0, 8.0, 8.0, 9.0, 2.0, 9.0, 6.0, 3.0, 3.0, 1.0, 3.0, 9.0, 3.0, 9.0, 1.0, 10.0, 3.0, 8.0, 10.0, 10.0, 6.0, 6.0, 3.0]
global b_y = 10
global p = [0.977, 0.155, 0.339, 0.732, 0.998, 0.378, 0.343, 0.581, 0.711, 0.217, 0.303, 0.805, 0.675, 0.705, 0.43, 0.898, 0.355, 0.794, 0.423, 0.941, 0.281, 0.832, 0.635, 0.855, 0.992, 0.427, 0.077, 0.392, 0.93, 0.885, 0.469, 0.395, 0.2, 0.491, 0.992, 0.562, 0.913, 0.575, 0.958, 0.859, 0.585, 0.124, 0.191, 0.155, 0.986, 0.066, 0.342, 0.421, 0.673, 0.928, 0.994, 0.797, 0.052, 0.947, 0.542, 0.54, 0.657, 0.243, 0.465, 0.989, 0.932, 0.675, 0.635, 0.271, 0.659, 0.958, 0.692, 0.4, 0.688, 0.683, 0.156, 0.256, 0.905, 0.346, 0.59, 0.036, 0.797, 0.741, 0.562, 0.653, 0.93, 0.409, 0.203, 0.081, 0.655, 0.386, 0.478, 0.176, 0.177, 0.034, 0.38, 0.493, 0.022, 0.643, 0.834, 0.71, 0.305, 0.744, 0.562, 0.289, 0.812, 0.601, 0.774, 0.566, 0.989, 0.108, 0.095, 0.147, 0.816, 0.33, 0.738, 0.113, 0.928, 0.669, 0.164, 0.664, 0.101, 0.778, 0.075, 0.349, 0.499, 0.462, 0.572, 0.157, 0.241, 0.142, 0.629, 0.4, 0.162, 0.368, 0.248, 0.47, 0.52, 0.535, 0.819, 0.321, 0.924, 0.941, 0.906, 0.948, 0.567, 0.36, 0.511, 0.21, 0.765, 0.884, 0.719, 0.923, 0.237, 0.979, 0.386, 0.707, 0.314, 0.326, 0.117, 0.082, 0.462, 0.708, 0.147, 0.961, 0.182, 0.548, 0.115, 0.163, 0.481, 0.452, 0.121, 0.136, 0.503, 0.812, 0.437, 0.923, 0.977, 0.065, 0.034, 0.534, 0.405, 0.468, 0.998, 0.369, 0.599, 0.844, 0.779, 0.485, 0.639, 0.38, 0.136, 0.841, 0.757, 0.162, 0.865, 0.567, 0.997, 0.421, 0.898, 0.965, 0.163, 0.009, 0.291, 0.743, 0.909, 0.13, 0.948, 0.168, 0.237, 0.104, 0.774, 0.409, 0.131, 0.764, 0.958, 0.357, 0.856, 0.347, 0.757, 0.128, 0.98, 0.119, 0.251, 0.07, 0.493, 0.709, 0.394, 0.281, 0.942, 0.72, 0.625, 0.008, 0.502, 0.58, 0.098, 0.47, 0.764, 0.117, 0.348, 0.454, 0.394, 0.89, 0.521, 0.096, 0.566, 0.365, 0.466, 0.676, 0.564, 0.747, 0.31, 0.472, 0.705, 0.869, 0.807, 0.864, 0.025, 0.926, 0.707, 0.463, 0.6, 0.254, 0.179, 0.575, 0.465, 0.597, 0.855, 0.927, 0.088, 0.911, 0.719, 0.788, 0.577, 0.88, 0.059, 0.518, 0.329, 0.665, 0.569, 0.115, 0.404, 0.199, 0.344, 0.005, 0.98, 0.699, 0.225, 0.212, 0.067, 0.808, 0.422, 0.465, 0.973, 0.355, 0.525, 0.25, 0.291, 0.028, 0.014, 0.497, 0.819, 0.866, 0.989, 0.362, 0.226, 0.017, 0.917, 0.362, 0.855, 0.025, 0.972, 0.398, 0.825, 0.516, 0.137, 0.873, 0.905, 0.374, 0.036, 0.838, 0.544, 0.824]
global q = [0.995, 0.747, 0.743, 0.996, 0.999, 0.477, 0.72, 0.924, 0.981, 0.328, 0.773, 0.865, 0.861, 0.97, 0.99, 0.935, 0.497, 0.934, 0.927, 0.981, 0.286, 0.945, 0.736, 0.903, 0.995, 0.571, 0.689, 0.83, 0.972, 0.915, 0.719, 0.967, 0.863, 0.53, 0.998, 0.899, 0.965, 0.835, 0.976, 0.926, 0.657, 0.917, 0.671, 0.186, 0.987, 0.148, 0.868, 0.912, 0.92, 0.995, 0.994, 0.994, 0.809, 0.981, 0.889, 0.969, 0.765, 0.763, 0.971, 0.99, 0.992, 0.864, 0.707, 0.704, 0.914, 0.987, 0.811, 0.488, 0.886, 0.747, 0.158, 0.712, 0.967, 0.605, 0.651, 0.911, 0.943, 0.886, 0.613, 0.822, 0.967, 0.565, 0.501, 0.767, 0.88, 0.647, 0.907, 0.745, 0.29, 0.667, 0.419, 0.985, 0.715, 0.695, 0.936, 0.787, 0.675, 0.91, 0.688, 0.344, 0.994, 0.932, 0.985, 0.631, 0.989, 0.225, 0.366, 0.149, 0.886, 0.394, 0.856, 0.62, 0.929, 0.922, 0.438, 0.786, 0.963, 0.987, 0.703, 0.777, 0.625, 0.534, 0.81, 0.315, 0.672, 0.244, 0.978, 0.641, 0.427, 0.755, 0.859, 0.614, 0.65, 0.677, 0.999, 0.98, 0.981, 0.961, 0.956, 0.948, 0.789, 0.804, 0.81, 0.263, 0.876, 0.998, 0.807, 0.978, 0.972, 0.984, 0.777, 0.753, 0.419, 0.92, 0.967, 0.908, 0.783, 0.924, 0.265, 0.993, 0.816, 0.944, 0.197, 0.167, 0.742, 0.933, 0.763, 0.205, 0.98, 0.941, 0.891, 0.984, 0.977, 0.194, 0.312, 0.895, 0.927, 0.544, 0.999, 0.685, 0.797, 0.847, 0.873, 0.899, 0.804, 0.467, 0.375, 0.91, 0.975, 0.634, 0.978, 0.681, 0.997, 0.814, 0.926, 0.978, 0.672, 0.482, 0.989, 0.743, 0.99, 0.912, 0.949, 0.595, 0.657, 0.849, 0.863, 0.897, 0.891, 0.929, 0.958, 0.935, 0.932, 0.354, 0.811, 0.481, 0.997, 0.691, 0.563, 0.285, 0.593, 0.82, 0.478, 0.703, 0.987, 0.74, 0.731, 0.393, 0.591, 0.626, 0.127, 0.887, 0.827, 0.415, 0.601, 0.649, 0.689, 0.911, 0.909, 0.166, 0.973, 0.465, 0.746, 0.828, 0.803, 0.948, 0.678, 0.496, 0.708, 0.947, 0.991, 0.931, 0.317, 0.931, 0.819, 0.534, 0.638, 0.355, 0.879, 0.724, 0.526, 0.947, 0.92, 0.951, 0.884, 0.952, 0.932, 0.931, 0.754, 0.881, 0.097, 0.926, 0.907, 0.743, 0.889, 0.575, 0.531, 0.677, 0.903, 0.059, 0.992, 0.723, 0.582, 0.965, 0.884, 0.945, 0.881, 0.794, 0.996, 0.438, 0.653, 0.744, 0.807, 0.808, 0.17, 0.865, 0.925, 0.887, 0.994, 0.692, 0.668, 0.747, 0.994, 0.657, 0.889, 0.094, 0.977, 0.565, 0.83, 0.912, 0.629, 0.904, 0.942, 0.903, 0.956, 0.983, 0.746, 0.89]
global origin = 1
global destination = 60