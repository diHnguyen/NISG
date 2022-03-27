global arcs = [1 14; 1 20; 1 24; 1 27; 1 31; 2 14; 2 45; 3 7; 3 11; 3 20; 3 21; 3 25; 3 36; 3 37; 3 40; 3 42; 3 44; 3 46; 3 48; 4 2; 4 6; 4 22; 4 41; 5 12; 5 24; 5 30; 5 35; 5 41; 5 42; 5 48; 5 49; 6 37; 7 8; 7 19; 8 26; 8 33; 8 38; 8 46; 9 4; 9 12; 9 21; 9 25; 9 28; 9 33; 9 35; 9 36; 9 41; 10 21; 10 23; 10 24; 10 26; 11 14; 11 28; 11 30; 11 35; 11 36; 11 38; 11 39; 11 46; 12 8; 12 13; 12 15; 12 24; 12 28; 12 32; 13 6; 13 12; 13 16; 13 21; 13 35; 13 43; 14 4; 14 6; 14 8; 14 13; 14 15; 14 22; 14 30; 14 41; 14 42; 14 43; 14 47; 15 4; 15 8; 15 12; 15 18; 15 21; 15 23; 15 26; 15 38; 16 9; 16 13; 16 29; 16 32; 16 46; 16 47; 16 49; 17 4; 17 6; 17 16; 17 20; 17 25; 17 31; 18 7; 18 10; 18 16; 18 27; 18 39; 19 2; 19 5; 19 12; 19 38; 19 47; 20 9; 20 18; 20 46; 21 30; 21 31; 21 32; 21 37; 21 41; 22 11; 22 27; 22 31; 22 36; 22 46; 23 3; 23 5; 23 7; 23 15; 23 19; 23 24; 23 31; 23 37; 24 3; 24 8; 24 15; 24 33; 24 38; 24 41; 24 46; 24 47; 25 10; 25 20; 25 22; 25 31; 25 33; 25 41; 26 3; 26 39; 26 46; 27 26; 27 30; 27 32; 27 35; 27 44; 28 7; 28 18; 28 31; 28 32; 29 5; 29 28; 29 37; 29 45; 29 49; 30 7; 30 9; 30 32; 30 44; 30 46; 31 8; 31 15; 31 17; 31 19; 31 30; 31 36; 31 40; 31 46; 31 48; 32 7; 32 17; 32 50; 33 34; 34 12; 34 13; 34 16; 34 25; 34 31; 34 42; 34 47; 34 48; 35 4; 35 10; 35 21; 35 34; 36 21; 36 25; 36 40; 36 47; 37 8; 37 9; 37 11; 37 24; 37 49; 38 14; 38 21; 38 23; 38 43; 38 49; 39 4; 39 5; 39 19; 39 33; 39 36; 39 50; 40 39; 41 5; 41 6; 41 18; 41 30; 41 32; 42 12; 42 22; 42 45; 43 4; 43 7; 43 21; 43 34; 43 38; 44 31; 44 45; 45 4; 45 17; 45 27; 45 28; 46 2; 46 4; 46 27; 46 38; 46 45; 47 19; 47 38; 48 13; 48 37; 48 50; 49 15; 49 17; 49 27]
global d_x = [2.0, 7.0, 2.0, 7.0, 4.0, 1.0, 6.0, 10.0, 8.0, 6.0, 5.0, 1.0, 2.0, 5.0, 1.0, 5.0, 8.0, 9.0, 2.0, 4.0, 2.0, 9.0, 2.0, 1.0, 9.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 1.0, 1.0, 8.0, 4.0, 1.0, 3.0, 3.0, 6.0, 5.0, 9.0, 2.0, 4.0, 8.0, 5.0, 5.0, 9.0, 9.0, 1.0, 1.0, 6.0, 6.0, 6.0, 9.0, 10.0, 3.0, 5.0, 2.0, 9.0, 2.0, 9.0, 5.0, 5.0, 1.0, 7.0, 4.0, 1.0, 6.0, 2.0, 1.0, 8.0, 2.0, 5.0, 8.0, 10.0, 9.0, 1.0, 9.0, 6.0, 4.0, 9.0, 8.0, 10.0, 6.0, 4.0, 8.0, 3.0, 4.0, 10.0, 4.0, 3.0, 2.0, 6.0, 3.0, 10.0, 6.0, 7.0, 4.0, 3.0, 8.0, 9.0, 4.0, 10.0, 5.0, 1.0, 6.0, 10.0, 6.0, 3.0, 10.0, 10.0, 3.0, 10.0, 6.0, 5.0, 8.0, 4.0, 3.0, 1.0, 5.0, 7.0, 10.0, 2.0, 9.0, 4.0, 4.0, 2.0, 7.0, 1.0, 4.0, 10.0, 1.0, 2.0, 8.0, 7.0, 1.0, 5.0, 1.0, 9.0, 10.0, 10.0, 2.0, 9.0, 9.0, 5.0, 3.0, 8.0, 1.0, 3.0, 5.0, 7.0, 1.0, 2.0, 7.0, 4.0, 3.0, 1.0, 8.0, 5.0, 3.0, 5.0, 10.0, 1.0, 4.0, 10.0, 7.0, 6.0, 1.0, 7.0, 6.0, 2.0, 4.0, 4.0, 4.0, 10.0, 2.0, 7.0, 7.0, 8.0, 10.0, 4.0, 2.0, 4.0, 10.0, 3.0, 3.0, 3.0, 9.0, 6.0, 10.0, 3.0, 3.0, 10.0, 8.0, 3.0, 10.0, 5.0, 4.0, 1.0, 3.0, 5.0, 9.0, 10.0, 4.0, 4.0, 9.0, 4.0, 1.0, 9.0, 8.0, 7.0, 1.0, 3.0, 7.0, 10.0, 7.0, 3.0, 3.0, 5.0, 3.0, 6.0, 1.0, 8.0, 6.0, 1.0, 2.0, 9.0, 3.0, 4.0, 6.0, 10.0, 9.0, 5.0, 3.0, 4.0, 6.0, 1.0, 1.0, 4.0, 4.0, 5.0, 8.0, 10.0, 3.0, 5.0, 8.0, 2.0, 4.0]
global b_x = 5
global d_y = [7.0, 7.0, 4.0, 9.0, 10.0, 4.0, 3.0, 10.0, 10.0, 9.0, 1.0, 5.0, 3.0, 2.0, 4.0, 1.0, 3.0, 3.0, 10.0, 2.0, 3.0, 1.0, 3.0, 4.0, 8.0, 7.0, 2.0, 5.0, 2.0, 6.0, 10.0, 5.0, 5.0, 8.0, 1.0, 1.0, 7.0, 6.0, 5.0, 3.0, 9.0, 5.0, 2.0, 6.0, 10.0, 8.0, 1.0, 1.0, 3.0, 9.0, 9.0, 3.0, 6.0, 1.0, 1.0, 5.0, 8.0, 2.0, 9.0, 4.0, 8.0, 5.0, 9.0, 6.0, 6.0, 6.0, 10.0, 3.0, 7.0, 3.0, 10.0, 4.0, 6.0, 6.0, 3.0, 9.0, 1.0, 3.0, 6.0, 9.0, 7.0, 8.0, 8.0, 2.0, 7.0, 8.0, 2.0, 10.0, 8.0, 7.0, 9.0, 6.0, 4.0, 9.0, 4.0, 7.0, 9.0, 7.0, 1.0, 7.0, 1.0, 1.0, 8.0, 8.0, 2.0, 9.0, 1.0, 9.0, 10.0, 3.0, 3.0, 10.0, 2.0, 10.0, 9.0, 5.0, 5.0, 9.0, 9.0, 2.0, 7.0, 4.0, 5.0, 8.0, 1.0, 1.0, 3.0, 9.0, 2.0, 1.0, 8.0, 7.0, 6.0, 2.0, 7.0, 10.0, 4.0, 3.0, 8.0, 2.0, 10.0, 3.0, 1.0, 9.0, 8.0, 5.0, 3.0, 7.0, 2.0, 10.0, 1.0, 6.0, 3.0, 4.0, 2.0, 7.0, 8.0, 2.0, 6.0, 6.0, 4.0, 4.0, 3.0, 2.0, 8.0, 1.0, 1.0, 5.0, 7.0, 1.0, 10.0, 6.0, 4.0, 10.0, 10.0, 6.0, 8.0, 5.0, 1.0, 1.0, 8.0, 4.0, 2.0, 6.0, 6.0, 8.0, 2.0, 5.0, 4.0, 9.0, 8.0, 7.0, 1.0, 1.0, 3.0, 2.0, 7.0, 5.0, 8.0, 10.0, 10.0, 10.0, 5.0, 6.0, 5.0, 8.0, 3.0, 4.0, 2.0, 4.0, 9.0, 10.0, 3.0, 3.0, 3.0, 10.0, 10.0, 3.0, 4.0, 9.0, 9.0, 6.0, 8.0, 3.0, 2.0, 9.0, 8.0, 9.0, 3.0, 9.0, 3.0, 1.0, 3.0, 4.0, 7.0, 1.0, 2.0, 4.0, 7.0, 6.0, 6.0, 8.0, 7.0, 7.0, 9.0, 8.0, 5.0, 2.0]
global b_y = 10
global p = [0.216, 0.832, 0.323, 0.912, 0.802, 0.841, 0.717, 0.328, 0.298, 0.156, 0.154, 0.314, 0.817, 0.34, 0.313, 0.608, 0.481, 0.557, 0.629, 0.358, 0.628, 0.148, 0.352, 0.029, 0.449, 0.754, 0.426, 0.299, 0.017, 0.189, 0.111, 0.18, 0.982, 0.305, 0.354, 0.166, 0.406, 0.939, 0.082, 0.701, 0.963, 0.904, 0.334, 0.313, 0.654, 0.499, 0.804, 0.207, 0.848, 0.721, 0.467, 0.585, 0.002, 0.698, 0.619, 0.809, 0.035, 0.706, 0.937, 0.189, 0.277, 0.469, 0.942, 0.515, 0.655, 0.623, 0.163, 0.623, 0.876, 0.536, 0.361, 0.182, 0.468, 0.043, 0.541, 0.807, 0.62, 0.554, 0.214, 0.963, 0.977, 0.244, 0.468, 0.197, 0.707, 0.332, 0.372, 0.569, 0.148, 0.29, 0.197, 0.972, 0.467, 0.769, 0.27, 0.97, 0.313, 0.481, 0.836, 0.112, 0.681, 0.025, 0.222, 0.92, 0.648, 0.666, 0.591, 0.31, 0.162, 0.185, 0.437, 0.777, 0.537, 0.921, 0.094, 0.488, 0.019, 0.492, 0.171, 0.994, 0.62, 0.525, 0.587, 0.46, 0.315, 0.388, 0.114, 0.933, 0.018, 0.272, 0.394, 0.789, 0.178, 0.486, 0.983, 0.63, 0.133, 0.388, 0.746, 0.415, 0.122, 0.54, 0.118, 0.828, 0.237, 0.448, 0.273, 0.804, 0.747, 0.386, 0.391, 0.833, 0.257, 0.661, 0.594, 0.159, 0.545, 0.543, 0.845, 0.828, 0.247, 0.153, 0.959, 0.919, 0.387, 0.434, 0.384, 0.022, 0.941, 0.704, 0.099, 0.909, 0.597, 0.693, 0.516, 0.229, 0.786, 0.615, 0.361, 0.254, 0.038, 0.332, 0.19, 0.611, 0.679, 0.81, 0.412, 0.76, 0.89, 0.017, 0.326, 0.358, 0.613, 0.838, 0.689, 0.545, 0.912, 0.537, 0.31, 0.372, 0.678, 0.213, 0.906, 0.412, 0.936, 0.033, 0.021, 0.128, 0.171, 0.966, 0.014, 0.522, 0.04, 0.584, 0.021, 0.563, 0.5, 0.099, 0.74, 0.066, 0.312, 0.062, 0.82, 0.915, 0.158, 0.253, 0.989, 0.194, 0.375, 0.303, 0.613, 0.16, 0.863, 0.592, 0.371, 0.367, 0.98, 0.143, 0.028, 0.487, 0.672, 0.578, 0.824, 0.049, 0.62, 0.838, 0.56, 0.668]
global q = [0.735, 0.944, 0.902, 0.969, 0.97, 0.961, 0.773, 0.692, 0.324, 0.69, 0.285, 0.497, 0.848, 0.657, 0.816, 0.979, 0.982, 0.966, 0.846, 0.535, 0.991, 0.695, 0.935, 0.42, 0.654, 0.854, 0.735, 0.71, 0.616, 0.819, 0.934, 0.278, 0.994, 0.796, 0.375, 0.801, 0.724, 0.985, 0.277, 0.986, 0.991, 0.92, 0.923, 0.387, 0.844, 0.643, 0.847, 0.222, 0.868, 0.753, 0.605, 0.609, 0.618, 0.861, 0.749, 0.944, 0.3, 0.787, 0.98, 0.335, 0.731, 0.725, 0.984, 0.955, 0.691, 0.836, 0.881, 0.801, 0.894, 0.977, 0.775, 0.347, 0.812, 0.34, 0.81, 0.969, 0.743, 0.998, 0.753, 0.992, 0.993, 0.536, 0.707, 0.658, 0.721, 0.489, 0.546, 0.665, 0.943, 0.733, 0.26, 0.997, 0.494, 0.818, 0.9, 0.984, 0.814, 0.669, 0.919, 0.146, 0.759, 0.564, 0.36, 0.986, 0.874, 0.827, 0.962, 0.623, 0.419, 0.865, 0.441, 0.787, 0.604, 0.993, 0.442, 0.682, 0.959, 0.582, 0.469, 0.999, 0.757, 0.61, 0.685, 0.984, 0.63, 0.455, 0.216, 0.994, 0.971, 0.41, 0.912, 0.893, 0.335, 0.799, 0.986, 0.77, 0.559, 0.667, 0.971, 0.885, 0.299, 0.658, 0.288, 0.987, 0.267, 0.6, 0.724, 0.901, 0.873, 0.908, 0.838, 0.897, 0.648, 0.779, 0.613, 0.895, 0.87, 0.695, 0.877, 0.946, 0.493, 0.602, 0.993, 0.974, 0.901, 0.634, 0.752, 0.461, 0.963, 0.969, 0.563, 0.945, 0.674, 0.801, 0.903, 0.586, 0.985, 0.882, 0.957, 0.487, 0.678, 0.403, 0.758, 0.777, 0.914, 0.84, 0.722, 0.831, 0.987, 0.395, 0.389, 0.898, 0.828, 0.973, 0.882, 0.66, 0.963, 0.709, 0.855, 0.629, 0.748, 0.44, 0.932, 0.513, 0.958, 0.242, 0.664, 0.697, 0.78, 0.996, 0.107, 0.612, 0.558, 0.828, 0.515, 0.972, 0.543, 0.848, 0.779, 0.44, 0.884, 0.264, 0.873, 0.953, 0.467, 0.686, 0.995, 0.364, 0.515, 0.378, 0.875, 0.791, 0.932, 0.685, 0.829, 0.499, 0.984, 0.438, 0.284, 0.861, 0.841, 0.845, 0.987, 0.579, 0.862, 0.905, 0.939, 0.999]
global origin = 1
global destination = 50