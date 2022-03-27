global arcs = [1 3; 1 9; 1 10; 1 12; 1 13; 1 17; 2 4; 2 5; 2 9; 2 10; 2 12; 2 15; 2 18; 2 20; 3 9; 3 15; 3 19; 4 2; 4 3; 4 7; 4 14; 4 16; 4 17; 4 18; 5 2; 5 8; 5 11; 5 15; 5 16; 5 19; 5 20; 6 4; 6 10; 6 11; 6 15; 6 16; 6 18; 6 20; 7 3; 7 10; 7 12; 7 14; 7 19; 7 20; 8 3; 8 10; 8 12; 8 13; 8 15; 9 3; 9 11; 9 13; 9 17; 9 20; 10 5; 10 6; 10 11; 10 14; 10 15; 10 18; 10 20; 11 2; 11 3; 11 10; 11 13; 11 18; 12 2; 12 3; 12 8; 13 7; 13 10; 13 12; 13 19; 14 2; 14 5; 14 6; 14 7; 14 15; 14 19; 15 3; 15 6; 15 10; 15 13; 15 16; 15 20; 16 3; 16 19; 17 2; 17 3; 17 16; 18 2; 18 3; 18 7; 18 10; 18 11; 18 13; 18 15; 18 17; 19 4; 19 7; 19 8; 19 11; 19 13; 19 15; 19 17; 19 18; 19 20]
global d_x = [1.0, 1.0, 10.0, 10.0, 7.0, 10.0, 2.0, 2.0, 3.0, 6.0, 5.0, 2.0, 2.0, 8.0, 6.0, 9.0, 5.0, 3.0, 6.0, 3.0, 7.0, 8.0, 4.0, 6.0, 10.0, 3.0, 10.0, 9.0, 10.0, 1.0, 2.0, 8.0, 4.0, 1.0, 9.0, 1.0, 3.0, 6.0, 10.0, 7.0, 9.0, 1.0, 9.0, 5.0, 7.0, 3.0, 9.0, 1.0, 4.0, 5.0, 1.0, 4.0, 4.0, 2.0, 4.0, 1.0, 4.0, 9.0, 10.0, 1.0, 10.0, 1.0, 7.0, 3.0, 9.0, 7.0, 3.0, 6.0, 2.0, 1.0, 10.0, 7.0, 6.0, 5.0, 7.0, 2.0, 7.0, 9.0, 5.0, 3.0, 8.0, 2.0, 9.0, 8.0, 4.0, 6.0, 10.0, 7.0, 5.0, 4.0, 8.0, 4.0, 5.0, 6.0, 1.0, 3.0, 3.0, 1.0, 4.0, 9.0, 6.0, 1.0, 9.0, 1.0, 4.0, 7.0, 7.0]
global b_x = 5
global d_y = [3.0, 1.0, 1.0, 2.0, 3.0, 8.0, 4.0, 1.0, 10.0, 3.0, 4.0, 2.0, 6.0, 6.0, 3.0, 3.0, 2.0, 2.0, 4.0, 5.0, 7.0, 6.0, 8.0, 8.0, 7.0, 3.0, 7.0, 10.0, 8.0, 1.0, 10.0, 10.0, 4.0, 1.0, 1.0, 8.0, 4.0, 8.0, 5.0, 9.0, 9.0, 3.0, 4.0, 7.0, 4.0, 7.0, 6.0, 1.0, 1.0, 6.0, 9.0, 7.0, 6.0, 1.0, 3.0, 1.0, 10.0, 1.0, 1.0, 1.0, 5.0, 1.0, 8.0, 8.0, 7.0, 4.0, 9.0, 7.0, 8.0, 7.0, 8.0, 6.0, 1.0, 7.0, 9.0, 7.0, 10.0, 8.0, 1.0, 6.0, 2.0, 7.0, 6.0, 9.0, 9.0, 2.0, 6.0, 10.0, 3.0, 10.0, 8.0, 7.0, 8.0, 2.0, 6.0, 1.0, 4.0, 3.0, 10.0, 5.0, 4.0, 10.0, 4.0, 3.0, 10.0, 8.0, 1.0]
global b_y = 10
global p = [0.269, 0.098, 0.054, 0.205, 0.512, 0.121, 0.852, 0.281, 0.783, 0.247, 0.392, 0.757, 0.189, 0.443, 0.591, 0.113, 0.542, 0.696, 0.869, 0.847, 0.235, 0.015, 0.911, 0.414, 0.358, 0.413, 0.776, 0.75, 0.892, 0.929, 0.081, 0.024, 0.717, 0.751, 0.467, 0.896, 0.28, 0.438, 0.182, 0.631, 0.545, 0.589, 0.557, 0.094, 0.359, 0.456, 0.588, 0.162, 0.892, 0.958, 0.175, 0.772, 0.428, 0.49, 0.467, 0.307, 0.727, 0.829, 0.289, 0.458, 0.765, 0.301, 0.024, 0.453, 0.995, 0.405, 0.896, 0.306, 0.231, 0.242, 0.875, 0.57, 0.165, 0.013, 0.703, 0.171, 0.894, 0.46, 0.073, 0.332, 0.68, 0.78, 0.983, 0.129, 0.151, 0.632, 0.57, 0.493, 0.513, 0.433, 0.319, 0.413, 0.371, 0.282, 0.505, 0.357, 0.667, 0.137, 0.063, 0.339, 0.927, 0.618, 0.591, 0.433, 0.798, 0.582, 0.345]
global q = [0.594, 0.108, 0.733, 0.409, 0.914, 0.753, 0.873, 0.922, 0.959, 0.931, 0.516, 0.905, 0.767, 0.747, 0.605, 0.645, 0.901, 0.783, 0.911, 0.853, 0.355, 0.134, 0.995, 0.871, 0.684, 0.851, 0.819, 0.914, 0.922, 0.992, 0.796, 0.183, 0.844, 0.892, 0.83, 0.92, 0.617, 0.971, 0.925, 0.757, 0.56, 0.809, 0.892, 0.573, 0.448, 0.913, 0.975, 0.227, 0.991, 0.968, 0.193, 0.943, 0.745, 0.817, 0.515, 0.842, 0.956, 0.848, 0.84, 0.524, 0.846, 0.678, 0.501, 0.691, 0.997, 0.699, 0.896, 0.495, 0.685, 0.312, 0.929, 0.644, 0.328, 0.864, 0.953, 0.924, 0.957, 0.656, 0.639, 0.556, 0.975, 0.939, 0.994, 0.523, 0.235, 0.843, 0.762, 0.607, 0.528, 0.708, 0.595, 0.571, 0.639, 0.61, 0.589, 0.438, 0.982, 0.228, 0.946, 0.915, 0.936, 0.752, 0.643, 0.697, 0.918, 0.6, 0.775]
global origin = 1
global destination = 20