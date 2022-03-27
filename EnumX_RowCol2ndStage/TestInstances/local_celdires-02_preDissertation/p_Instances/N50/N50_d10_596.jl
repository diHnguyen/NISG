global arcs = [1 29; 1 32; 1 39; 1 40; 2 3; 2 9; 2 23; 2 24; 2 38; 3 9; 3 17; 3 25; 4 20; 4 22; 5 3; 5 37; 5 39; 6 12; 6 17; 6 27; 6 29; 6 31; 7 2; 7 8; 7 9; 7 12; 7 21; 7 30; 7 32; 7 39; 7 44; 7 50; 8 3; 8 4; 8 5; 8 39; 8 47; 9 3; 9 13; 9 15; 9 45; 10 6; 10 7; 10 8; 10 14; 10 18; 10 21; 10 29; 10 35; 11 2; 11 17; 11 40; 11 44; 11 45; 12 30; 12 44; 12 46; 13 11; 13 21; 13 24; 14 11; 14 37; 14 38; 14 43; 15 36; 15 46; 15 47; 16 44; 16 45; 17 12; 17 18; 17 22; 17 24; 17 26; 17 43; 18 22; 18 32; 18 36; 18 47; 19 2; 19 6; 19 9; 19 14; 19 17; 19 24; 19 36; 19 46; 19 47; 20 2; 20 15; 20 18; 20 19; 20 26; 20 27; 21 11; 21 28; 21 29; 21 30; 21 34; 21 38; 21 42; 22 8; 22 9; 22 15; 22 20; 22 31; 22 32; 22 33; 22 41; 23 27; 23 35; 23 39; 23 45; 24 4; 24 14; 25 3; 25 12; 25 13; 25 19; 25 20; 25 23; 25 30; 25 39; 25 48; 26 13; 26 29; 26 38; 26 40; 26 42; 27 26; 27 28; 27 29; 27 36; 27 37; 27 38; 27 41; 28 17; 28 34; 28 36; 29 22; 29 41; 29 44; 29 47; 30 14; 30 17; 31 36; 31 40; 32 6; 32 7; 32 11; 32 28; 32 45; 33 12; 33 15; 33 29; 33 40; 34 2; 34 4; 34 25; 34 38; 35 8; 35 13; 35 36; 35 38; 35 43; 36 7; 36 15; 36 19; 36 32; 36 49; 36 50; 37 34; 37 38; 38 8; 38 13; 38 14; 38 15; 38 17; 38 18; 38 27; 38 49; 38 50; 39 19; 39 20; 39 22; 39 27; 39 30; 40 30; 40 33; 40 47; 41 16; 42 7; 42 15; 42 21; 42 34; 43 10; 43 20; 43 23; 43 44; 44 6; 44 9; 44 10; 44 28; 44 31; 44 34; 44 43; 44 47; 45 3; 45 23; 45 31; 45 36; 45 38; 46 20; 46 22; 46 31; 46 47; 47 8; 47 11; 47 14; 47 16; 47 17; 47 34; 47 36; 47 37; 47 46; 48 23; 48 28; 48 32; 48 36; 48 40; 48 42; 48 46; 48 49; 49 7]
global d_x = [7.0, 2.0, 2.0, 10.0, 6.0, 2.0, 3.0, 2.0, 8.0, 6.0, 7.0, 4.0, 1.0, 1.0, 4.0, 8.0, 3.0, 3.0, 6.0, 8.0, 4.0, 5.0, 1.0, 6.0, 7.0, 8.0, 1.0, 2.0, 8.0, 5.0, 7.0, 4.0, 3.0, 2.0, 6.0, 10.0, 5.0, 2.0, 4.0, 1.0, 9.0, 5.0, 5.0, 1.0, 10.0, 6.0, 7.0, 6.0, 3.0, 8.0, 9.0, 9.0, 7.0, 7.0, 8.0, 10.0, 7.0, 7.0, 4.0, 6.0, 5.0, 3.0, 10.0, 8.0, 5.0, 8.0, 8.0, 4.0, 10.0, 3.0, 2.0, 6.0, 8.0, 10.0, 9.0, 8.0, 1.0, 9.0, 8.0, 5.0, 5.0, 1.0, 1.0, 3.0, 9.0, 4.0, 2.0, 10.0, 5.0, 7.0, 5.0, 9.0, 1.0, 6.0, 4.0, 1.0, 2.0, 2.0, 2.0, 10.0, 4.0, 7.0, 5.0, 8.0, 7.0, 3.0, 7.0, 3.0, 9.0, 3.0, 10.0, 9.0, 7.0, 4.0, 4.0, 2.0, 6.0, 9.0, 7.0, 2.0, 5.0, 9.0, 8.0, 2.0, 10.0, 2.0, 7.0, 6.0, 9.0, 2.0, 3.0, 2.0, 6.0, 4.0, 2.0, 10.0, 9.0, 5.0, 5.0, 3.0, 6.0, 5.0, 6.0, 5.0, 10.0, 1.0, 5.0, 1.0, 8.0, 8.0, 4.0, 9.0, 8.0, 2.0, 9.0, 10.0, 3.0, 2.0, 10.0, 6.0, 5.0, 9.0, 1.0, 1.0, 9.0, 4.0, 4.0, 4.0, 4.0, 6.0, 7.0, 8.0, 8.0, 5.0, 6.0, 7.0, 9.0, 3.0, 5.0, 7.0, 1.0, 5.0, 2.0, 4.0, 2.0, 2.0, 5.0, 5.0, 4.0, 8.0, 2.0, 5.0, 9.0, 2.0, 1.0, 5.0, 10.0, 3.0, 6.0, 5.0, 8.0, 6.0, 3.0, 3.0, 3.0, 6.0, 8.0, 7.0, 7.0, 10.0, 8.0, 9.0, 10.0, 5.0, 8.0, 9.0, 7.0, 9.0, 8.0, 1.0, 1.0, 10.0, 3.0, 5.0, 7.0, 6.0, 3.0, 7.0, 9.0, 8.0, 6.0, 8.0, 9.0, 10.0]
global b_x = 5
global d_y = [10.0, 2.0, 7.0, 9.0, 4.0, 9.0, 4.0, 9.0, 7.0, 9.0, 10.0, 8.0, 10.0, 7.0, 6.0, 9.0, 4.0, 1.0, 9.0, 6.0, 1.0, 8.0, 6.0, 10.0, 10.0, 8.0, 9.0, 2.0, 7.0, 5.0, 10.0, 5.0, 1.0, 3.0, 9.0, 8.0, 7.0, 10.0, 3.0, 9.0, 3.0, 8.0, 2.0, 10.0, 6.0, 6.0, 3.0, 6.0, 10.0, 6.0, 2.0, 1.0, 2.0, 5.0, 8.0, 1.0, 4.0, 8.0, 9.0, 5.0, 5.0, 9.0, 8.0, 4.0, 6.0, 3.0, 4.0, 8.0, 6.0, 6.0, 1.0, 4.0, 1.0, 10.0, 10.0, 8.0, 5.0, 9.0, 2.0, 3.0, 1.0, 10.0, 8.0, 1.0, 1.0, 2.0, 7.0, 3.0, 7.0, 4.0, 9.0, 9.0, 6.0, 4.0, 4.0, 4.0, 6.0, 2.0, 2.0, 4.0, 4.0, 8.0, 2.0, 3.0, 4.0, 10.0, 10.0, 10.0, 10.0, 4.0, 3.0, 10.0, 9.0, 9.0, 3.0, 6.0, 3.0, 2.0, 4.0, 6.0, 6.0, 7.0, 1.0, 8.0, 7.0, 6.0, 5.0, 9.0, 5.0, 4.0, 8.0, 6.0, 6.0, 5.0, 4.0, 6.0, 1.0, 10.0, 3.0, 4.0, 6.0, 3.0, 3.0, 7.0, 6.0, 5.0, 8.0, 4.0, 4.0, 6.0, 6.0, 7.0, 7.0, 6.0, 2.0, 6.0, 6.0, 4.0, 9.0, 1.0, 4.0, 10.0, 5.0, 5.0, 3.0, 7.0, 7.0, 1.0, 10.0, 5.0, 6.0, 4.0, 9.0, 7.0, 8.0, 9.0, 6.0, 7.0, 1.0, 7.0, 1.0, 3.0, 6.0, 5.0, 10.0, 6.0, 8.0, 10.0, 1.0, 8.0, 4.0, 7.0, 6.0, 6.0, 5.0, 6.0, 1.0, 1.0, 6.0, 6.0, 5.0, 5.0, 1.0, 8.0, 2.0, 10.0, 2.0, 10.0, 8.0, 9.0, 6.0, 8.0, 6.0, 8.0, 1.0, 9.0, 10.0, 1.0, 2.0, 5.0, 2.0, 4.0, 8.0, 10.0, 10.0, 5.0, 5.0, 2.0, 9.0, 4.0, 5.0, 1.0, 5.0, 5.0]
global b_y = 10
global p = [0.27, 0.663, 0.309, 0.583, 0.375, 0.106, 0.209, 0.198, 0.004, 0.011, 0.75, 0.315, 0.014, 0.531, 0.148, 0.109, 0.763, 0.465, 0.424, 0.427, 0.37, 0.834, 0.275, 0.429, 0.656, 0.512, 0.742, 0.345, 0.387, 0.258, 0.283, 0.953, 0.434, 0.288, 0.712, 0.67, 0.399, 0.425, 0.236, 0.907, 0.158, 0.797, 0.986, 0.563, 0.641, 0.685, 0.138, 0.466, 0.974, 0.12, 0.006, 0.817, 0.283, 0.827, 0.792, 0.678, 0.374, 0.958, 0.283, 0.109, 0.477, 0.869, 0.594, 0.093, 0.519, 0.535, 0.096, 0.546, 0.264, 0.478, 0.291, 0.567, 0.242, 0.704, 0.905, 0.262, 0.574, 0.026, 0.769, 0.995, 0.225, 0.835, 0.94, 0.471, 0.9, 0.639, 0.619, 0.337, 0.583, 0.346, 0.786, 0.875, 0.866, 0.505, 0.555, 0.514, 0.281, 0.512, 0.756, 0.637, 0.628, 0.344, 0.908, 0.512, 0.517, 0.995, 0.72, 0.713, 0.902, 0.594, 0.788, 0.674, 0.864, 0.01, 0.195, 0.987, 0.773, 0.294, 0.309, 0.677, 0.808, 0.373, 0.977, 0.511, 0.882, 0.156, 0.645, 0.224, 0.737, 0.875, 0.696, 0.027, 0.871, 0.944, 0.888, 0.928, 0.149, 0.425, 0.562, 0.76, 0.684, 0.602, 0.832, 0.123, 0.308, 0.98, 0.874, 0.573, 0.03, 0.71, 0.044, 0.901, 0.706, 0.336, 0.569, 0.112, 0.77, 0.74, 0.279, 0.22, 0.556, 0.586, 0.541, 0.629, 0.757, 0.049, 0.147, 0.231, 0.039, 0.679, 0.935, 0.486, 0.356, 0.76, 0.243, 0.496, 0.219, 0.725, 0.875, 0.016, 0.003, 0.218, 0.67, 0.798, 0.916, 0.91, 0.471, 0.739, 0.227, 0.747, 0.774, 0.258, 0.533, 0.23, 0.089, 0.769, 0.084, 0.944, 0.563, 0.925, 0.912, 0.912, 0.365, 0.546, 0.166, 0.257, 0.194, 0.25, 0.729, 0.225, 0.427, 0.234, 0.252, 0.703, 0.864, 0.905, 0.51, 0.037, 0.051, 0.444, 0.993, 0.058, 0.102, 0.285, 0.401, 0.465, 0.304, 0.977, 0.171, 0.95, 0.903, 0.575, 0.514, 0.505]
global q = [0.82, 0.919, 0.367, 0.824, 0.655, 0.442, 0.821, 0.865, 0.583, 0.685, 0.834, 0.472, 0.4, 0.932, 0.55, 0.747, 0.964, 0.582, 0.469, 0.958, 0.506, 0.936, 0.321, 0.445, 0.885, 0.836, 0.747, 0.616, 0.547, 0.944, 0.446, 0.957, 0.88, 0.57, 0.825, 0.796, 0.554, 0.481, 0.436, 0.954, 0.512, 0.825, 0.997, 0.663, 0.657, 0.985, 0.319, 0.898, 0.997, 0.384, 0.052, 0.824, 0.667, 0.999, 0.805, 0.694, 0.576, 0.975, 0.497, 0.677, 0.589, 0.96, 0.848, 0.757, 0.656, 0.996, 0.203, 0.914, 0.574, 0.627, 0.606, 0.904, 0.421, 0.774, 0.997, 0.971, 0.575, 0.449, 0.815, 0.999, 0.855, 0.962, 0.951, 0.946, 0.901, 0.749, 0.85, 0.881, 0.966, 0.773, 0.832, 0.984, 0.915, 0.996, 0.616, 0.74, 0.572, 0.517, 0.887, 0.948, 0.742, 0.79, 0.994, 0.57, 0.712, 0.995, 0.8, 0.831, 0.917, 0.975, 0.969, 0.891, 0.967, 0.912, 0.241, 0.997, 0.974, 0.991, 0.823, 0.715, 0.906, 0.431, 0.984, 0.914, 0.967, 0.874, 0.809, 0.461, 0.752, 0.947, 0.96, 0.905, 0.902, 0.948, 0.897, 0.995, 0.699, 0.438, 0.834, 0.925, 0.786, 0.652, 0.917, 0.987, 0.334, 0.995, 0.989, 0.7, 0.084, 0.955, 0.162, 0.978, 0.797, 0.848, 0.993, 0.596, 0.887, 0.802, 0.951, 0.802, 0.724, 0.855, 0.655, 0.77, 0.767, 0.389, 0.673, 0.848, 0.203, 0.736, 0.971, 0.576, 0.943, 0.832, 0.651, 0.713, 0.766, 0.935, 0.918, 0.393, 0.579, 0.731, 0.679, 0.827, 0.928, 0.949, 0.86, 0.806, 0.886, 0.758, 0.934, 0.841, 0.994, 0.766, 0.72, 0.874, 0.915, 0.945, 0.832, 0.981, 0.982, 0.92, 0.982, 0.942, 0.947, 0.342, 0.875, 0.599, 0.986, 0.372, 0.875, 0.789, 0.676, 0.768, 0.964, 0.985, 0.906, 0.229, 0.256, 0.543, 0.997, 0.572, 0.669, 0.497, 0.482, 0.553, 0.866, 0.997, 0.885, 0.966, 0.912, 0.885, 0.89, 0.862]
global origin = 1
global destination = 50