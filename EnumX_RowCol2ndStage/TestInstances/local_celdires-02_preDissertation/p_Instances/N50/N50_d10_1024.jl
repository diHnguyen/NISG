global arcs = [1 2; 1 16; 1 28; 1 37; 1 41; 2 28; 2 31; 2 40; 3 2; 3 8; 3 14; 3 17; 3 21; 3 23; 3 26; 3 39; 3 41; 3 43; 3 44; 3 49; 4 3; 4 14; 4 15; 4 39; 5 13; 5 45; 5 48; 6 16; 6 20; 6 26; 6 27; 6 30; 6 40; 6 44; 7 2; 7 23; 7 24; 7 34; 7 40; 7 44; 7 47; 7 49; 8 2; 8 15; 8 18; 8 31; 8 46; 8 49; 9 14; 9 37; 10 4; 10 16; 10 26; 10 43; 10 45; 10 46; 11 7; 11 15; 11 17; 11 39; 11 49; 12 5; 12 11; 12 25; 12 35; 12 47; 13 4; 13 7; 13 20; 13 30; 13 35; 13 44; 14 7; 14 11; 14 18; 14 30; 14 32; 14 46; 14 47; 14 50; 15 5; 15 12; 15 25; 15 35; 15 45; 15 48; 16 5; 16 6; 16 17; 16 27; 17 5; 17 11; 17 13; 17 22; 17 29; 17 35; 17 39; 17 46; 18 3; 18 5; 18 14; 18 16; 18 20; 18 23; 18 25; 18 36; 18 45; 18 47; 19 11; 19 15; 19 27; 19 30; 19 36; 19 40; 20 15; 20 18; 20 21; 20 30; 20 49; 21 7; 21 12; 21 19; 21 20; 21 30; 21 37; 21 47; 22 27; 22 29; 22 37; 22 43; 22 47; 23 6; 23 7; 23 15; 24 11; 24 12; 24 35; 24 47; 25 3; 25 7; 25 40; 26 9; 26 16; 26 19; 26 24; 26 27; 27 14; 27 42; 27 45; 27 48; 28 6; 28 17; 28 18; 28 27; 28 34; 28 41; 28 46; 29 6; 29 9; 29 13; 29 23; 29 24; 29 44; 30 20; 30 34; 30 37; 30 39; 31 4; 31 9; 31 11; 31 13; 31 27; 31 36; 32 4; 32 5; 32 10; 32 23; 32 24; 32 39; 32 47; 32 48; 33 3; 33 10; 33 29; 33 35; 33 38; 34 20; 34 30; 34 31; 35 2; 35 13; 35 22; 36 9; 36 20; 36 23; 37 19; 37 23; 37 27; 37 30; 37 36; 38 13; 38 36; 38 37; 39 35; 39 36; 40 3; 40 6; 40 11; 40 23; 40 33; 40 35; 40 47; 41 3; 41 12; 41 25; 41 27; 41 34; 41 35; 42 2; 42 47; 43 3; 43 12; 43 15; 43 36; 44 23; 44 27; 44 29; 44 35; 45 8; 45 10; 45 22; 45 32; 45 48; 45 49; 46 10; 46 19; 46 33; 46 37; 46 45; 47 11; 47 13; 47 28; 47 29; 47 33; 47 38; 48 9; 48 15; 48 16; 48 19; 48 25; 48 28; 48 36; 49 46; 49 50]
global d_x = [2.0, 5.0, 10.0, 10.0, 5.0, 7.0, 5.0, 1.0, 3.0, 6.0, 10.0, 5.0, 5.0, 7.0, 3.0, 9.0, 7.0, 4.0, 6.0, 9.0, 6.0, 8.0, 10.0, 5.0, 10.0, 1.0, 2.0, 5.0, 7.0, 10.0, 8.0, 3.0, 2.0, 3.0, 4.0, 3.0, 5.0, 9.0, 5.0, 1.0, 6.0, 5.0, 2.0, 8.0, 5.0, 5.0, 5.0, 6.0, 3.0, 3.0, 6.0, 6.0, 10.0, 1.0, 6.0, 2.0, 9.0, 2.0, 3.0, 3.0, 3.0, 6.0, 9.0, 6.0, 1.0, 1.0, 10.0, 7.0, 3.0, 4.0, 9.0, 1.0, 7.0, 10.0, 9.0, 10.0, 3.0, 2.0, 6.0, 4.0, 5.0, 8.0, 3.0, 5.0, 4.0, 3.0, 6.0, 1.0, 3.0, 2.0, 10.0, 10.0, 6.0, 2.0, 9.0, 10.0, 9.0, 1.0, 6.0, 6.0, 7.0, 10.0, 3.0, 8.0, 3.0, 5.0, 6.0, 2.0, 10.0, 7.0, 7.0, 1.0, 10.0, 7.0, 6.0, 8.0, 2.0, 9.0, 6.0, 6.0, 4.0, 1.0, 2.0, 2.0, 8.0, 10.0, 10.0, 9.0, 1.0, 4.0, 9.0, 6.0, 9.0, 8.0, 8.0, 10.0, 1.0, 10.0, 3.0, 5.0, 10.0, 10.0, 6.0, 9.0, 4.0, 6.0, 4.0, 8.0, 7.0, 9.0, 10.0, 8.0, 6.0, 1.0, 4.0, 7.0, 9.0, 9.0, 2.0, 2.0, 5.0, 6.0, 6.0, 5.0, 9.0, 2.0, 3.0, 7.0, 3.0, 3.0, 6.0, 7.0, 5.0, 3.0, 2.0, 7.0, 5.0, 9.0, 1.0, 1.0, 8.0, 5.0, 3.0, 10.0, 5.0, 5.0, 10.0, 4.0, 2.0, 7.0, 3.0, 1.0, 6.0, 7.0, 8.0, 10.0, 3.0, 3.0, 8.0, 6.0, 3.0, 5.0, 5.0, 3.0, 1.0, 9.0, 1.0, 3.0, 2.0, 2.0, 5.0, 1.0, 4.0, 10.0, 4.0, 5.0, 9.0, 10.0, 7.0, 8.0, 9.0, 10.0, 6.0, 5.0, 8.0, 8.0, 1.0, 7.0, 3.0, 10.0, 2.0, 8.0, 7.0, 8.0, 7.0, 4.0, 2.0, 7.0, 3.0, 2.0, 8.0, 5.0, 8.0, 3.0, 10.0, 1.0, 7.0, 9.0, 8.0, 8.0, 7.0, 3.0, 8.0, 6.0]
global b_x = 5
global d_y = [6.0, 10.0, 1.0, 1.0, 7.0, 7.0, 3.0, 8.0, 2.0, 9.0, 3.0, 1.0, 7.0, 7.0, 7.0, 6.0, 2.0, 1.0, 7.0, 9.0, 1.0, 4.0, 9.0, 9.0, 8.0, 3.0, 7.0, 7.0, 3.0, 10.0, 10.0, 2.0, 6.0, 2.0, 10.0, 2.0, 10.0, 3.0, 10.0, 9.0, 4.0, 1.0, 6.0, 3.0, 10.0, 4.0, 1.0, 8.0, 1.0, 7.0, 3.0, 2.0, 4.0, 4.0, 7.0, 5.0, 1.0, 2.0, 5.0, 2.0, 8.0, 9.0, 6.0, 8.0, 3.0, 6.0, 4.0, 10.0, 4.0, 8.0, 4.0, 9.0, 8.0, 4.0, 6.0, 5.0, 7.0, 8.0, 5.0, 6.0, 2.0, 10.0, 6.0, 2.0, 10.0, 1.0, 8.0, 3.0, 6.0, 10.0, 7.0, 4.0, 5.0, 5.0, 3.0, 8.0, 5.0, 7.0, 8.0, 7.0, 2.0, 6.0, 5.0, 9.0, 5.0, 8.0, 10.0, 5.0, 9.0, 7.0, 6.0, 1.0, 4.0, 8.0, 6.0, 9.0, 6.0, 4.0, 2.0, 7.0, 2.0, 2.0, 7.0, 5.0, 10.0, 1.0, 4.0, 3.0, 10.0, 9.0, 4.0, 7.0, 6.0, 1.0, 9.0, 6.0, 5.0, 5.0, 9.0, 3.0, 4.0, 1.0, 3.0, 10.0, 5.0, 9.0, 5.0, 7.0, 3.0, 8.0, 2.0, 4.0, 2.0, 8.0, 6.0, 5.0, 5.0, 10.0, 4.0, 2.0, 8.0, 2.0, 5.0, 1.0, 5.0, 10.0, 8.0, 7.0, 3.0, 9.0, 6.0, 4.0, 7.0, 8.0, 7.0, 6.0, 7.0, 9.0, 9.0, 1.0, 9.0, 7.0, 6.0, 2.0, 5.0, 4.0, 7.0, 10.0, 5.0, 7.0, 2.0, 1.0, 10.0, 1.0, 9.0, 7.0, 3.0, 10.0, 8.0, 2.0, 4.0, 7.0, 5.0, 9.0, 10.0, 5.0, 9.0, 9.0, 8.0, 8.0, 2.0, 5.0, 3.0, 3.0, 2.0, 8.0, 7.0, 3.0, 2.0, 9.0, 7.0, 3.0, 7.0, 4.0, 8.0, 8.0, 10.0, 5.0, 8.0, 10.0, 5.0, 4.0, 10.0, 10.0, 6.0, 5.0, 9.0, 9.0, 8.0, 3.0, 10.0, 10.0, 2.0, 2.0, 1.0, 8.0, 1.0, 1.0, 9.0, 5.0, 4.0, 3.0, 7.0, 2.0]
global b_y = 10
global p = [0.936, 0.533, 0.881, 0.231, 0.166, 0.081, 0.824, 0.759, 0.586, 0.041, 0.367, 0.768, 0.309, 0.327, 0.824, 0.236, 0.371, 0.37, 0.4, 0.535, 0.691, 0.83, 0.24, 0.195, 0.973, 0.962, 0.946, 0.721, 0.856, 0.214, 0.596, 0.722, 0.821, 0.972, 0.321, 0.241, 0.517, 0.508, 0.455, 0.093, 0.615, 0.691, 0.434, 0.861, 0.16, 0.564, 0.534, 0.048, 0.379, 0.072, 0.658, 0.28, 0.252, 0.092, 0.713, 0.692, 0.074, 0.908, 0.445, 0.246, 0.627, 0.19, 0.703, 0.078, 0.979, 0.677, 0.786, 0.07, 0.413, 0.326, 0.283, 0.066, 0.821, 0.607, 0.242, 0.161, 0.365, 0.084, 0.241, 0.311, 0.233, 0.475, 0.585, 0.159, 0.863, 0.901, 0.897, 0.04, 0.37, 0.197, 0.547, 0.995, 0.996, 0.042, 0.976, 0.421, 0.93, 0.086, 0.929, 0.283, 0.417, 0.001, 0.063, 0.73, 0.601, 0.469, 0.738, 0.516, 0.577, 0.859, 0.493, 0.78, 0.633, 0.033, 0.995, 0.455, 0.125, 0.241, 0.423, 0.322, 0.724, 0.236, 0.088, 0.205, 0.538, 0.277, 0.154, 0.391, 0.281, 0.908, 0.151, 0.157, 0.451, 0.733, 0.854, 0.363, 0.246, 0.107, 0.077, 0.686, 0.747, 0.023, 0.693, 0.447, 0.564, 0.571, 0.825, 0.212, 0.666, 0.498, 0.285, 0.611, 0.596, 0.172, 0.521, 0.434, 0.756, 0.402, 0.122, 0.889, 0.871, 0.957, 0.319, 0.485, 0.274, 0.379, 0.224, 0.553, 0.729, 0.98, 0.959, 0.499, 0.601, 0.462, 0.48, 0.697, 0.787, 0.149, 0.222, 0.125, 0.455, 0.883, 0.333, 0.607, 0.133, 0.055, 0.804, 0.66, 0.89, 0.869, 0.265, 0.482, 0.689, 0.736, 0.93, 0.776, 0.871, 0.799, 0.511, 0.885, 0.941, 0.879, 0.82, 0.161, 0.135, 0.631, 0.457, 0.765, 0.013, 0.02, 0.557, 0.049, 0.239, 0.117, 0.784, 0.107, 0.774, 0.392, 0.298, 0.893, 0.084, 0.522, 0.354, 0.566, 0.749, 0.053, 0.315, 0.694, 0.799, 0.885, 0.988, 0.524, 0.505, 0.495, 0.653, 0.092, 0.394, 0.674, 0.686, 0.691, 0.773, 0.336, 0.083, 0.391, 0.046, 0.677, 0.762, 0.794, 0.552, 0.878, 0.921, 0.532, 0.292, 0.051]
global q = [0.974, 0.708, 0.908, 0.455, 0.74, 0.138, 0.886, 0.895, 0.907, 0.964, 0.385, 0.967, 0.436, 0.464, 0.826, 0.345, 0.431, 0.541, 0.431, 0.847, 0.726, 0.96, 0.906, 0.497, 0.995, 0.983, 0.951, 0.774, 0.912, 0.375, 0.787, 0.816, 0.873, 0.974, 0.971, 0.468, 0.868, 0.908, 0.545, 0.727, 0.975, 0.7, 0.755, 0.988, 0.353, 0.946, 0.562, 0.856, 0.768, 0.396, 0.864, 0.626, 0.285, 0.262, 0.716, 0.733, 0.163, 0.951, 0.544, 0.714, 0.975, 0.595, 0.984, 0.677, 0.991, 0.857, 0.955, 0.353, 0.453, 0.592, 0.811, 0.673, 0.874, 0.921, 0.684, 0.203, 0.716, 0.978, 0.496, 0.994, 0.691, 0.67, 0.974, 0.621, 0.966, 0.992, 0.982, 0.986, 0.721, 0.505, 0.805, 0.995, 0.998, 0.617, 0.981, 0.631, 0.995, 0.184, 0.943, 0.873, 0.932, 0.003, 0.459, 0.735, 0.715, 0.72, 0.803, 0.747, 0.775, 0.887, 0.649, 0.801, 0.71, 0.404, 0.999, 0.666, 0.421, 0.292, 0.888, 0.54, 0.739, 0.918, 0.961, 0.847, 0.854, 0.462, 0.343, 0.857, 0.344, 0.96, 0.277, 0.252, 0.673, 0.965, 0.979, 0.99, 0.296, 0.655, 0.57, 0.887, 0.978, 0.817, 0.849, 0.82, 0.734, 0.704, 0.839, 0.708, 0.948, 0.697, 0.789, 0.636, 0.609, 0.35, 0.885, 0.923, 0.764, 0.609, 0.687, 0.944, 0.874, 0.981, 0.831, 0.507, 0.474, 0.805, 0.933, 0.72, 0.845, 0.988, 0.978, 0.917, 0.893, 0.945, 0.782, 0.787, 0.966, 0.698, 0.991, 0.714, 0.987, 0.903, 0.34, 0.935, 0.317, 0.14, 0.829, 0.816, 0.909, 0.977, 0.343, 0.882, 0.86, 0.908, 0.98, 0.99, 0.983, 0.912, 0.81, 0.981, 0.962, 0.983, 0.971, 0.844, 0.741, 0.917, 0.563, 0.969, 0.026, 0.764, 0.671, 0.834, 0.482, 0.828, 0.853, 0.849, 0.936, 0.718, 0.714, 0.928, 0.932, 0.748, 0.638, 0.695, 0.788, 0.719, 0.74, 0.832, 0.824, 0.906, 0.99, 0.871, 0.996, 0.866, 0.722, 0.184, 0.83, 0.782, 0.844, 0.985, 0.907, 0.827, 0.413, 0.644, 0.928, 0.95, 0.85, 0.906, 0.736, 0.984, 0.929, 0.99, 0.614, 0.186]
global origin = 1
global destination = 50