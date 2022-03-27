global arcs = [1 2; 1 4; 1 17; 1 19; 1 29; 1 45; 2 7; 2 25; 2 40; 2 47; 3 2; 3 4; 3 30; 3 42; 4 2; 4 7; 4 14; 4 26; 4 49; 5 9; 5 36; 6 5; 6 17; 6 44; 7 6; 7 8; 7 18; 7 20; 7 44; 8 3; 8 7; 8 19; 8 37; 8 47; 9 6; 9 7; 9 18; 9 20; 9 25; 9 36; 9 40; 9 45; 9 47; 9 49; 10 28; 10 45; 11 7; 11 15; 11 24; 11 39; 11 40; 12 2; 12 16; 12 17; 12 33; 13 2; 13 7; 13 9; 13 21; 13 46; 13 49; 14 43; 15 16; 15 21; 15 48; 16 2; 16 12; 16 18; 16 23; 16 27; 17 19; 17 21; 17 29; 17 42; 17 43; 17 46; 18 12; 18 28; 18 38; 18 43; 19 20; 19 23; 19 34; 19 37; 20 23; 20 39; 20 41; 21 9; 21 19; 21 24; 21 26; 21 30; 21 44; 22 7; 22 25; 22 26; 22 38; 22 43; 22 47; 23 6; 23 26; 23 28; 23 46; 24 16; 24 18; 24 38; 24 40; 24 43; 24 48; 25 20; 25 26; 25 27; 25 29; 25 31; 25 36; 25 38; 25 42; 26 7; 26 9; 26 15; 26 32; 26 33; 27 11; 27 25; 27 26; 27 39; 27 41; 28 4; 28 6; 28 11; 28 16; 28 17; 28 27; 28 44; 29 13; 29 14; 29 15; 29 18; 29 28; 29 37; 29 47; 29 50; 30 2; 30 12; 30 28; 30 31; 30 35; 30 37; 30 46; 31 3; 31 9; 31 18; 31 34; 31 47; 32 2; 32 11; 32 18; 32 26; 32 38; 32 39; 32 48; 33 32; 33 40; 33 46; 34 16; 34 19; 34 21; 34 30; 35 12; 35 17; 35 23; 35 37; 36 4; 36 9; 36 20; 36 23; 36 35; 36 38; 37 18; 37 36; 37 39; 38 12; 38 24; 38 35; 38 44; 38 45; 39 3; 39 21; 40 2; 40 6; 40 46; 41 7; 41 19; 41 22; 41 33; 41 34; 41 36; 41 44; 42 17; 43 7; 43 20; 43 22; 43 37; 43 44; 43 48; 44 7; 44 18; 44 19; 44 25; 44 27; 44 36; 44 47; 45 31; 45 44; 46 3; 46 25; 46 30; 46 39; 47 10; 47 11; 47 21; 47 39; 47 45; 48 6; 48 19; 48 27; 48 41; 49 9; 49 13; 49 29; 49 50]
global d_x = [7.0, 4.0, 8.0, 9.0, 6.0, 9.0, 3.0, 7.0, 10.0, 7.0, 3.0, 7.0, 7.0, 1.0, 10.0, 3.0, 7.0, 4.0, 1.0, 5.0, 5.0, 4.0, 3.0, 1.0, 5.0, 9.0, 5.0, 1.0, 8.0, 4.0, 8.0, 5.0, 8.0, 8.0, 3.0, 5.0, 8.0, 2.0, 3.0, 7.0, 2.0, 8.0, 10.0, 2.0, 8.0, 6.0, 8.0, 1.0, 10.0, 4.0, 9.0, 3.0, 2.0, 1.0, 4.0, 9.0, 4.0, 8.0, 5.0, 3.0, 7.0, 3.0, 2.0, 9.0, 7.0, 1.0, 6.0, 2.0, 3.0, 10.0, 3.0, 7.0, 7.0, 10.0, 7.0, 3.0, 6.0, 5.0, 8.0, 7.0, 4.0, 3.0, 6.0, 2.0, 1.0, 4.0, 4.0, 7.0, 6.0, 8.0, 1.0, 2.0, 6.0, 3.0, 9.0, 5.0, 1.0, 10.0, 8.0, 7.0, 4.0, 2.0, 10.0, 5.0, 9.0, 9.0, 2.0, 3.0, 4.0, 2.0, 5.0, 6.0, 10.0, 2.0, 8.0, 10.0, 5.0, 7.0, 6.0, 10.0, 6.0, 3.0, 8.0, 2.0, 2.0, 4.0, 1.0, 7.0, 10.0, 8.0, 1.0, 8.0, 4.0, 5.0, 5.0, 8.0, 9.0, 4.0, 10.0, 9.0, 10.0, 8.0, 8.0, 1.0, 4.0, 1.0, 4.0, 3.0, 1.0, 2.0, 10.0, 2.0, 1.0, 5.0, 1.0, 6.0, 10.0, 9.0, 8.0, 5.0, 3.0, 1.0, 3.0, 9.0, 10.0, 9.0, 1.0, 2.0, 2.0, 9.0, 6.0, 10.0, 4.0, 6.0, 2.0, 1.0, 9.0, 2.0, 3.0, 1.0, 8.0, 3.0, 7.0, 5.0, 10.0, 1.0, 4.0, 7.0, 6.0, 8.0, 6.0, 9.0, 6.0, 4.0, 7.0, 7.0, 8.0, 6.0, 4.0, 7.0, 1.0, 8.0, 7.0, 4.0, 5.0, 8.0, 6.0, 6.0, 6.0, 3.0, 5.0, 1.0, 3.0, 10.0, 4.0, 10.0, 1.0, 2.0, 6.0, 10.0, 8.0, 1.0, 3.0, 5.0, 6.0, 6.0, 8.0, 10.0, 4.0, 6.0, 8.0]
global b_x = 5
global d_y = [8.0, 6.0, 3.0, 4.0, 9.0, 3.0, 9.0, 6.0, 3.0, 6.0, 10.0, 3.0, 4.0, 8.0, 1.0, 2.0, 10.0, 2.0, 7.0, 4.0, 2.0, 4.0, 2.0, 9.0, 2.0, 5.0, 1.0, 10.0, 1.0, 3.0, 8.0, 7.0, 4.0, 1.0, 9.0, 3.0, 10.0, 3.0, 7.0, 4.0, 3.0, 5.0, 3.0, 6.0, 4.0, 5.0, 9.0, 10.0, 7.0, 2.0, 3.0, 4.0, 10.0, 1.0, 2.0, 1.0, 10.0, 9.0, 8.0, 10.0, 5.0, 4.0, 6.0, 7.0, 4.0, 6.0, 1.0, 5.0, 3.0, 3.0, 7.0, 1.0, 10.0, 6.0, 8.0, 6.0, 6.0, 6.0, 2.0, 5.0, 4.0, 4.0, 7.0, 10.0, 7.0, 5.0, 8.0, 4.0, 3.0, 8.0, 5.0, 1.0, 4.0, 5.0, 5.0, 2.0, 6.0, 1.0, 5.0, 1.0, 3.0, 1.0, 9.0, 3.0, 6.0, 9.0, 8.0, 2.0, 1.0, 6.0, 9.0, 10.0, 8.0, 10.0, 6.0, 5.0, 5.0, 5.0, 8.0, 9.0, 9.0, 5.0, 8.0, 9.0, 10.0, 6.0, 3.0, 3.0, 6.0, 9.0, 9.0, 3.0, 3.0, 4.0, 1.0, 1.0, 7.0, 3.0, 10.0, 2.0, 5.0, 2.0, 6.0, 5.0, 5.0, 6.0, 1.0, 5.0, 8.0, 9.0, 4.0, 2.0, 2.0, 10.0, 9.0, 3.0, 4.0, 10.0, 7.0, 7.0, 4.0, 1.0, 8.0, 7.0, 7.0, 4.0, 5.0, 6.0, 5.0, 7.0, 8.0, 7.0, 1.0, 1.0, 9.0, 2.0, 8.0, 10.0, 4.0, 1.0, 9.0, 1.0, 1.0, 6.0, 5.0, 9.0, 9.0, 6.0, 6.0, 6.0, 8.0, 9.0, 4.0, 2.0, 9.0, 6.0, 2.0, 8.0, 2.0, 5.0, 8.0, 9.0, 7.0, 9.0, 4.0, 9.0, 10.0, 5.0, 10.0, 6.0, 5.0, 5.0, 3.0, 1.0, 9.0, 9.0, 3.0, 1.0, 6.0, 7.0, 2.0, 1.0, 1.0, 5.0, 5.0, 1.0, 1.0, 9.0, 7.0, 1.0, 3.0]
global b_y = 10
global p = [0.024, 0.707, 0.759, 0.146, 0.956, 0.293, 0.619, 0.32, 0.712, 0.356, 0.159, 0.659, 0.417, 0.073, 0.254, 0.427, 0.464, 0.557, 0.55, 0.041, 0.353, 0.826, 0.139, 0.6, 0.561, 0.6, 0.126, 0.224, 0.349, 0.258, 0.451, 0.939, 0.757, 0.406, 0.598, 0.158, 0.674, 0.08, 0.047, 0.509, 0.52, 0.891, 0.03, 0.906, 0.961, 0.54, 0.898, 0.759, 0.876, 0.387, 0.689, 0.815, 0.938, 0.261, 0.909, 0.105, 0.344, 0.819, 0.671, 0.483, 0.557, 0.786, 0.617, 0.766, 0.638, 0.336, 0.021, 0.949, 0.942, 0.713, 0.831, 0.892, 0.224, 0.757, 0.386, 0.886, 0.312, 0.568, 0.96, 0.971, 0.334, 0.068, 0.555, 0.409, 0.77, 0.484, 0.193, 0.092, 0.179, 0.537, 0.351, 0.604, 0.416, 0.64, 0.17, 0.576, 0.105, 0.809, 0.031, 0.005, 0.254, 0.396, 0.628, 0.667, 0.385, 0.941, 0.313, 0.3, 0.503, 0.263, 0.651, 0.725, 0.297, 0.575, 0.495, 0.789, 0.621, 0.475, 0.142, 0.659, 0.587, 0.76, 0.636, 0.662, 0.023, 0.469, 0.345, 0.33, 0.792, 0.711, 0.784, 0.547, 0.231, 0.345, 0.201, 0.13, 0.298, 0.844, 0.54, 0.671, 0.252, 0.492, 0.962, 0.39, 0.967, 0.028, 0.465, 0.508, 0.027, 0.034, 0.672, 0.272, 0.554, 0.918, 0.23, 0.938, 0.188, 0.245, 0.217, 0.114, 0.711, 0.523, 0.268, 0.558, 0.359, 0.945, 0.021, 0.085, 0.665, 0.383, 0.361, 0.996, 0.661, 0.851, 0.018, 0.617, 0.785, 0.287, 0.842, 0.679, 0.097, 0.358, 0.01, 0.345, 0.273, 0.496, 0.262, 0.41, 0.879, 0.426, 0.839, 0.099, 0.124, 0.924, 0.925, 0.269, 0.125, 0.224, 0.36, 0.875, 0.952, 0.486, 0.095, 0.528, 0.643, 0.429, 0.907, 0.141, 0.986, 0.808, 0.662, 0.178, 0.902, 0.494, 0.265, 0.416, 0.266, 0.667, 0.588, 0.551, 0.111, 0.93, 0.792, 0.537, 0.284, 0.55, 0.403, 0.842, 0.83, 0.928, 0.486]
global q = [0.912, 0.839, 0.809, 0.7, 0.96, 0.661, 0.812, 0.882, 0.717, 0.586, 0.312, 0.785, 0.873, 0.513, 0.368, 0.867, 0.995, 0.734, 0.964, 0.945, 0.637, 0.929, 0.37, 0.733, 0.59, 0.628, 0.464, 0.278, 0.836, 0.905, 0.835, 0.985, 0.917, 0.581, 0.977, 0.253, 0.819, 0.277, 0.358, 0.699, 0.74, 0.958, 0.296, 0.958, 0.984, 0.556, 0.92, 0.941, 0.938, 0.846, 0.921, 0.844, 0.982, 0.491, 0.958, 0.484, 0.494, 0.96, 0.78, 0.753, 0.592, 0.975, 0.856, 0.774, 0.778, 0.934, 0.248, 0.983, 0.998, 0.741, 0.864, 0.902, 0.478, 0.799, 0.932, 0.998, 0.453, 0.647, 0.975, 0.98, 0.392, 0.526, 0.597, 0.476, 0.968, 0.867, 0.452, 0.134, 0.885, 0.944, 0.398, 0.964, 0.654, 0.868, 0.903, 0.787, 0.71, 0.999, 0.523, 0.723, 0.768, 0.488, 0.781, 0.87, 0.966, 0.994, 0.994, 0.76, 0.785, 0.709, 0.708, 0.89, 0.467, 0.852, 0.526, 0.817, 0.7, 0.965, 0.574, 0.84, 0.926, 0.799, 0.926, 0.664, 0.096, 0.682, 0.579, 0.712, 0.884, 0.75, 0.922, 0.929, 0.614, 0.887, 0.701, 0.935, 0.598, 0.906, 0.854, 0.851, 0.802, 0.74, 0.976, 0.507, 0.975, 0.54, 0.968, 0.934, 0.998, 0.043, 0.768, 0.913, 0.837, 0.948, 0.319, 0.98, 0.328, 0.659, 0.281, 0.66, 0.848, 0.644, 0.913, 0.907, 0.81, 0.957, 0.189, 0.613, 0.815, 0.435, 0.619, 0.996, 0.893, 0.941, 0.995, 0.988, 0.834, 0.686, 0.976, 0.887, 0.557, 0.732, 0.304, 0.799, 0.999, 0.779, 0.741, 0.442, 0.992, 0.982, 0.895, 0.133, 0.943, 0.96, 0.937, 0.473, 0.979, 0.361, 0.683, 0.896, 0.982, 0.855, 0.706, 0.696, 0.94, 0.829, 0.986, 0.904, 0.989, 0.994, 0.944, 0.41, 0.975, 0.851, 0.521, 0.654, 0.809, 0.909, 0.748, 0.936, 0.309, 0.943, 0.947, 0.724, 0.853, 0.892, 0.901, 0.991, 0.898, 0.961, 0.945]
global origin = 1
global destination = 50