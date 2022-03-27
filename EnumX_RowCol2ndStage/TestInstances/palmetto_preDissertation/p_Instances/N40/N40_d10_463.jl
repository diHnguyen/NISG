global arcs = [1 5; 1 9; 1 11; 1 34; 2 4; 2 7; 2 19; 2 30; 2 35; 2 36; 3 24; 3 29; 3 35; 3 38; 4 23; 5 6; 5 10; 5 11; 5 25; 5 28; 6 27; 6 28; 7 3; 7 8; 7 11; 7 20; 7 23; 7 25; 7 35; 7 38; 8 17; 8 19; 8 35; 9 13; 9 17; 9 26; 9 32; 10 2; 10 7; 10 8; 10 18; 10 30; 10 34; 10 35; 11 20; 12 30; 13 5; 13 6; 13 7; 13 9; 13 20; 13 22; 13 24; 13 32; 14 2; 14 12; 14 19; 14 23; 15 2; 15 13; 15 14; 15 21; 15 36; 15 38; 16 4; 16 17; 16 19; 16 33; 16 38; 17 3; 17 6; 18 17; 19 25; 19 29; 19 40; 20 4; 20 5; 20 22; 20 36; 20 37; 21 15; 21 16; 22 5; 22 15; 22 20; 22 30; 23 2; 23 29; 23 32; 23 35; 23 40; 24 5; 24 12; 24 27; 24 31; 25 3; 25 13; 25 23; 25 28; 25 36; 26 29; 26 30; 26 40; 27 14; 27 36; 28 10; 28 19; 28 30; 28 31; 29 10; 29 25; 30 5; 30 7; 30 9; 30 10; 30 24; 31 14; 31 28; 31 39; 32 7; 32 31; 32 36; 33 4; 33 5; 33 11; 33 23; 34 15; 34 25; 34 33; 34 35; 35 14; 35 26; 35 28; 35 30; 35 33; 36 10; 36 12; 36 28; 36 34; 37 5; 37 15; 37 19; 37 22; 37 25; 37 30; 37 40; 38 5; 38 8; 38 11; 38 13; 38 14; 38 15; 38 16; 38 31; 39 14; 39 21]
global d_x = [3.0, 2.0, 9.0, 8.0, 1.0, 2.0, 3.0, 8.0, 3.0, 5.0, 1.0, 1.0, 5.0, 8.0, 6.0, 9.0, 1.0, 3.0, 5.0, 10.0, 1.0, 2.0, 6.0, 8.0, 3.0, 6.0, 4.0, 1.0, 9.0, 2.0, 8.0, 6.0, 3.0, 3.0, 1.0, 10.0, 6.0, 4.0, 6.0, 6.0, 3.0, 5.0, 6.0, 6.0, 2.0, 2.0, 1.0, 2.0, 6.0, 6.0, 1.0, 8.0, 9.0, 8.0, 10.0, 6.0, 2.0, 9.0, 8.0, 7.0, 3.0, 7.0, 10.0, 10.0, 4.0, 3.0, 6.0, 4.0, 7.0, 8.0, 9.0, 4.0, 9.0, 1.0, 9.0, 10.0, 6.0, 8.0, 2.0, 2.0, 9.0, 8.0, 1.0, 4.0, 4.0, 2.0, 3.0, 1.0, 8.0, 6.0, 9.0, 4.0, 10.0, 8.0, 1.0, 4.0, 5.0, 5.0, 7.0, 6.0, 10.0, 7.0, 4.0, 5.0, 1.0, 5.0, 5.0, 5.0, 10.0, 6.0, 10.0, 7.0, 6.0, 2.0, 5.0, 7.0, 3.0, 6.0, 3.0, 2.0, 9.0, 6.0, 10.0, 2.0, 3.0, 1.0, 6.0, 9.0, 2.0, 8.0, 8.0, 4.0, 1.0, 3.0, 4.0, 8.0, 6.0, 5.0, 3.0, 4.0, 5.0, 9.0, 8.0, 7.0, 4.0, 5.0, 9.0, 5.0, 5.0, 2.0, 3.0, 1.0, 5.0, 2.0, 8.0, 8.0]
global b_x = 5
global d_y = [7.0, 7.0, 6.0, 5.0, 1.0, 6.0, 7.0, 10.0, 8.0, 1.0, 9.0, 2.0, 10.0, 5.0, 8.0, 5.0, 9.0, 8.0, 6.0, 3.0, 4.0, 1.0, 2.0, 1.0, 2.0, 5.0, 2.0, 9.0, 7.0, 6.0, 4.0, 3.0, 9.0, 1.0, 8.0, 3.0, 10.0, 5.0, 4.0, 1.0, 4.0, 5.0, 8.0, 8.0, 7.0, 3.0, 4.0, 8.0, 3.0, 3.0, 4.0, 8.0, 10.0, 9.0, 8.0, 2.0, 3.0, 7.0, 7.0, 10.0, 9.0, 7.0, 9.0, 3.0, 7.0, 4.0, 4.0, 10.0, 3.0, 1.0, 7.0, 5.0, 10.0, 6.0, 1.0, 5.0, 5.0, 3.0, 4.0, 1.0, 8.0, 4.0, 5.0, 5.0, 4.0, 1.0, 9.0, 9.0, 8.0, 7.0, 9.0, 9.0, 5.0, 6.0, 2.0, 6.0, 10.0, 9.0, 5.0, 10.0, 3.0, 4.0, 7.0, 4.0, 7.0, 8.0, 5.0, 4.0, 1.0, 9.0, 2.0, 7.0, 6.0, 9.0, 7.0, 1.0, 6.0, 10.0, 9.0, 4.0, 7.0, 5.0, 3.0, 10.0, 8.0, 1.0, 1.0, 5.0, 7.0, 2.0, 2.0, 6.0, 6.0, 3.0, 7.0, 10.0, 4.0, 1.0, 10.0, 9.0, 3.0, 6.0, 8.0, 4.0, 10.0, 10.0, 10.0, 1.0, 1.0, 4.0, 8.0, 7.0, 7.0, 8.0, 9.0, 5.0]
global b_y = 10
global p = [0.805, 0.92, 0.674, 0.61, 0.072, 0.602, 0.615, 0.324, 0.468, 0.892, 0.71, 0.715, 0.447, 0.875, 0.756, 0.625, 0.455, 0.878, 0.226, 0.682, 0.162, 0.895, 0.261, 0.127, 0.021, 0.767, 0.857, 0.593, 0.436, 0.408, 0.164, 0.533, 0.766, 0.257, 0.445, 0.568, 0.627, 0.682, 0.224, 0.99, 0.895, 0.527, 0.156, 0.145, 0.492, 0.332, 0.246, 0.444, 0.542, 0.73, 0.639, 0.123, 0.132, 0.92, 0.968, 0.491, 0.343, 0.23, 0.477, 0.073, 0.439, 0.403, 0.605, 0.267, 0.223, 0.396, 0.824, 0.419, 0.423, 0.72, 0.633, 0.658, 0.506, 0.551, 0.26, 0.014, 0.169, 0.842, 0.912, 0.523, 0.435, 0.584, 0.073, 0.256, 0.294, 0.037, 0.034, 0.098, 0.674, 0.252, 0.723, 0.349, 0.301, 0.065, 0.061, 0.264, 0.321, 0.706, 0.278, 0.195, 0.77, 0.739, 0.464, 0.188, 0.443, 0.939, 0.274, 0.103, 0.781, 0.999, 0.633, 0.801, 0.846, 0.185, 0.185, 0.997, 0.878, 0.834, 0.909, 0.374, 0.199, 0.889, 0.432, 0.336, 0.643, 0.289, 0.332, 0.233, 0.21, 0.373, 0.846, 0.349, 0.008, 0.31, 0.831, 0.523, 0.247, 0.315, 0.812, 0.267, 0.067, 0.873, 0.463, 0.596, 0.587, 0.146, 0.168, 0.582, 0.793, 0.916, 0.56, 0.722, 0.108, 0.195, 0.865, 0.468]
global q = [0.838, 0.959, 0.743, 0.81, 0.097, 0.946, 0.813, 0.99, 0.559, 0.976, 0.931, 0.817, 0.816, 0.923, 0.863, 0.948, 0.748, 0.944, 0.967, 0.98, 0.312, 0.951, 0.57, 0.196, 0.327, 0.959, 0.888, 0.873, 0.939, 0.549, 0.816, 0.637, 0.81, 0.457, 0.511, 0.704, 0.791, 0.866, 0.316, 0.991, 0.949, 0.819, 0.445, 0.476, 0.752, 0.641, 0.541, 0.703, 0.592, 0.771, 0.737, 0.459, 0.413, 0.922, 0.972, 0.562, 0.87, 0.306, 0.797, 0.444, 0.873, 0.659, 0.959, 0.428, 0.52, 0.416, 0.872, 0.76, 0.778, 0.899, 0.864, 0.876, 0.616, 0.626, 0.935, 0.071, 0.315, 0.991, 0.989, 0.853, 0.555, 0.793, 0.271, 0.483, 0.742, 0.908, 0.546, 0.561, 0.762, 0.345, 0.943, 0.657, 0.78, 0.368, 0.835, 0.827, 0.748, 0.876, 0.968, 0.868, 0.81, 0.879, 0.926, 0.231, 0.455, 0.978, 0.576, 0.819, 0.95, 0.999, 0.837, 0.891, 0.856, 0.966, 0.342, 0.999, 0.884, 0.948, 0.984, 0.833, 0.937, 0.952, 0.517, 0.667, 0.932, 0.779, 0.728, 0.336, 0.823, 0.393, 0.992, 0.745, 0.349, 0.57, 0.893, 0.792, 0.516, 0.724, 0.952, 0.954, 0.57, 0.895, 0.819, 0.814, 0.721, 0.456, 0.447, 0.991, 0.877, 0.975, 0.579, 0.807, 0.921, 0.304, 0.93, 0.71]
global origin = 1
global destination = 40