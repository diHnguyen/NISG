global arcs = [1 2; 1 3; 1 14; 1 25; 1 28; 1 32; 1 47; 2 6; 2 8; 2 11; 2 12; 2 15; 2 18; 2 29; 2 32; 2 37; 2 43; 3 21; 3 23; 3 47; 4 8; 4 12; 4 20; 4 22; 4 24; 4 31; 4 32; 4 48; 5 13; 5 27; 5 34; 5 38; 5 43; 5 49; 6 11; 6 20; 6 28; 6 30; 6 35; 6 42; 6 44; 7 9; 7 13; 7 20; 7 23; 7 25; 7 26; 7 30; 7 42; 7 45; 8 3; 8 15; 9 11; 9 15; 9 16; 9 22; 9 34; 9 37; 9 40; 10 6; 10 11; 10 21; 10 23; 10 27; 10 30; 10 32; 10 44; 11 6; 11 24; 11 32; 11 35; 12 33; 12 38; 12 45; 13 5; 13 9; 13 16; 13 43; 14 5; 14 13; 14 18; 14 29; 15 3; 15 4; 15 5; 15 12; 15 17; 15 24; 15 25; 16 2; 16 12; 16 18; 16 19; 16 21; 16 38; 16 39; 17 37; 17 43; 18 13; 18 25; 18 26; 18 33; 18 50; 19 36; 19 43; 19 46; 20 3; 20 18; 20 27; 20 36; 20 45; 20 46; 20 48; 21 13; 21 23; 21 24; 21 30; 22 20; 22 33; 22 43; 23 15; 23 35; 23 39; 23 43; 23 48; 23 50; 24 35; 24 36; 24 38; 24 39; 25 7; 25 8; 25 12; 25 13; 25 26; 25 33; 25 45; 26 14; 26 15; 26 34; 26 43; 27 5; 27 22; 27 48; 28 5; 28 49; 29 19; 29 22; 29 23; 29 24; 29 26; 29 27; 29 32; 29 47; 30 4; 30 8; 30 10; 30 17; 30 26; 30 28; 30 36; 31 9; 31 10; 31 26; 31 41; 31 44; 32 6; 32 9; 32 26; 32 30; 32 38; 32 40; 32 42; 32 44; 33 2; 33 22; 34 2; 34 10; 34 27; 34 35; 34 40; 35 11; 35 19; 35 30; 35 31; 35 37; 35 40; 36 22; 36 24; 36 30; 36 31; 36 33; 36 43; 36 44; 37 3; 37 10; 37 22; 37 25; 38 14; 38 48; 38 50; 39 8; 39 13; 39 17; 39 42; 39 47; 40 4; 40 41; 41 13; 41 15; 41 29; 41 38; 41 47; 42 2; 42 44; 43 2; 43 8; 43 21; 43 28; 43 30; 43 39; 44 12; 44 25; 45 4; 45 12; 45 15; 45 34; 45 36; 46 11; 46 27; 46 28; 46 32; 46 37; 46 39; 47 15; 47 20; 47 24; 47 29; 47 37; 47 40; 48 10; 48 15; 48 20; 48 31; 48 43; 48 50; 49 16; 49 23; 49 24; 49 30; 49 33; 49 34]
global d_x = [5.0, 8.0, 9.0, 9.0, 6.0, 9.0, 5.0, 6.0, 7.0, 10.0, 3.0, 9.0, 5.0, 9.0, 1.0, 9.0, 4.0, 10.0, 9.0, 5.0, 8.0, 7.0, 1.0, 9.0, 6.0, 10.0, 1.0, 10.0, 3.0, 1.0, 6.0, 10.0, 1.0, 1.0, 10.0, 3.0, 5.0, 9.0, 8.0, 8.0, 3.0, 3.0, 4.0, 2.0, 8.0, 6.0, 1.0, 8.0, 3.0, 8.0, 9.0, 9.0, 3.0, 4.0, 5.0, 10.0, 8.0, 6.0, 8.0, 10.0, 4.0, 8.0, 3.0, 10.0, 9.0, 8.0, 8.0, 2.0, 2.0, 5.0, 5.0, 9.0, 9.0, 3.0, 5.0, 10.0, 9.0, 4.0, 3.0, 10.0, 7.0, 9.0, 5.0, 1.0, 1.0, 5.0, 6.0, 1.0, 1.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 9.0, 2.0, 3.0, 2.0, 6.0, 3.0, 8.0, 10.0, 2.0, 9.0, 2.0, 2.0, 10.0, 8.0, 9.0, 6.0, 8.0, 5.0, 6.0, 7.0, 4.0, 9.0, 3.0, 10.0, 2.0, 5.0, 1.0, 1.0, 10.0, 8.0, 4.0, 2.0, 4.0, 5.0, 10.0, 8.0, 10.0, 10.0, 2.0, 7.0, 2.0, 6.0, 10.0, 7.0, 1.0, 4.0, 2.0, 7.0, 5.0, 9.0, 10.0, 5.0, 3.0, 9.0, 6.0, 8.0, 3.0, 10.0, 9.0, 3.0, 2.0, 6.0, 10.0, 6.0, 3.0, 5.0, 6.0, 10.0, 3.0, 10.0, 2.0, 7.0, 4.0, 8.0, 2.0, 10.0, 9.0, 8.0, 3.0, 5.0, 10.0, 3.0, 3.0, 5.0, 5.0, 6.0, 7.0, 2.0, 1.0, 9.0, 8.0, 3.0, 7.0, 6.0, 5.0, 1.0, 9.0, 10.0, 3.0, 4.0, 10.0, 4.0, 6.0, 4.0, 5.0, 6.0, 6.0, 2.0, 4.0, 1.0, 5.0, 7.0, 3.0, 4.0, 6.0, 2.0, 5.0, 7.0, 1.0, 9.0, 9.0, 5.0, 10.0, 5.0, 3.0, 8.0, 8.0, 1.0, 8.0, 8.0, 1.0, 7.0, 2.0, 4.0, 3.0, 10.0, 1.0, 5.0, 4.0, 5.0, 6.0, 9.0, 6.0, 10.0, 3.0, 1.0, 10.0, 3.0, 7.0, 4.0, 2.0, 8.0, 4.0, 2.0, 10.0, 6.0, 3.0]
global b_x = 5
global d_y = [6.0, 8.0, 8.0, 4.0, 4.0, 5.0, 5.0, 9.0, 1.0, 4.0, 2.0, 8.0, 7.0, 7.0, 10.0, 4.0, 2.0, 1.0, 4.0, 1.0, 6.0, 4.0, 5.0, 6.0, 8.0, 1.0, 8.0, 5.0, 9.0, 7.0, 10.0, 5.0, 9.0, 7.0, 9.0, 5.0, 6.0, 5.0, 3.0, 9.0, 7.0, 3.0, 2.0, 7.0, 4.0, 9.0, 10.0, 8.0, 6.0, 3.0, 1.0, 5.0, 4.0, 10.0, 5.0, 2.0, 4.0, 3.0, 6.0, 6.0, 9.0, 8.0, 10.0, 1.0, 6.0, 6.0, 6.0, 3.0, 9.0, 6.0, 3.0, 3.0, 6.0, 7.0, 4.0, 5.0, 9.0, 9.0, 10.0, 7.0, 7.0, 3.0, 2.0, 9.0, 10.0, 8.0, 4.0, 3.0, 9.0, 1.0, 6.0, 8.0, 7.0, 2.0, 10.0, 6.0, 7.0, 7.0, 6.0, 3.0, 8.0, 2.0, 7.0, 1.0, 1.0, 1.0, 9.0, 6.0, 5.0, 2.0, 7.0, 5.0, 1.0, 9.0, 1.0, 2.0, 9.0, 1.0, 1.0, 10.0, 6.0, 7.0, 1.0, 7.0, 7.0, 2.0, 7.0, 9.0, 2.0, 8.0, 8.0, 4.0, 7.0, 4.0, 1.0, 10.0, 8.0, 6.0, 10.0, 1.0, 2.0, 4.0, 8.0, 10.0, 5.0, 4.0, 4.0, 6.0, 8.0, 8.0, 5.0, 5.0, 4.0, 5.0, 10.0, 5.0, 10.0, 8.0, 1.0, 5.0, 4.0, 2.0, 10.0, 3.0, 8.0, 4.0, 5.0, 4.0, 9.0, 2.0, 6.0, 5.0, 6.0, 4.0, 2.0, 5.0, 3.0, 4.0, 5.0, 2.0, 6.0, 8.0, 3.0, 5.0, 8.0, 7.0, 9.0, 9.0, 5.0, 5.0, 7.0, 3.0, 8.0, 1.0, 2.0, 6.0, 4.0, 3.0, 2.0, 9.0, 7.0, 5.0, 3.0, 2.0, 10.0, 3.0, 3.0, 4.0, 3.0, 1.0, 7.0, 5.0, 6.0, 7.0, 1.0, 6.0, 9.0, 3.0, 2.0, 8.0, 10.0, 1.0, 10.0, 1.0, 2.0, 10.0, 10.0, 10.0, 7.0, 7.0, 2.0, 8.0, 5.0, 4.0, 7.0, 5.0, 7.0, 2.0, 5.0, 8.0, 5.0, 6.0, 7.0, 8.0, 7.0, 10.0, 2.0, 4.0, 5.0, 9.0, 8.0, 1.0]
global b_y = 10
global p = [0.138, 0.407, 0.579, 0.558, 0.396, 0.547, 0.09, 0.289, 0.038, 0.962, 0.219, 0.908, 0.73, 0.694, 0.236, 0.156, 0.054, 0.051, 0.854, 0.089, 0.239, 0.543, 0.005, 0.772, 0.98, 0.029, 0.466, 0.214, 0.608, 0.174, 0.379, 0.164, 0.753, 0.698, 0.193, 0.998, 0.073, 0.549, 0.675, 0.591, 0.221, 0.01, 0.317, 0.042, 0.225, 0.609, 0.313, 0.768, 0.923, 0.184, 0.547, 0.629, 0.904, 0.104, 0.907, 0.204, 0.233, 0.511, 0.426, 0.492, 0.209, 0.436, 0.166, 0.127, 0.903, 0.139, 0.897, 0.433, 0.603, 0.305, 0.049, 0.203, 0.043, 0.871, 0.203, 0.528, 0.103, 0.71, 0.198, 0.035, 0.889, 0.741, 0.206, 0.945, 0.409, 0.444, 0.728, 0.983, 0.517, 0.515, 0.387, 0.783, 0.169, 0.909, 0.405, 0.212, 0.917, 0.424, 0.606, 0.154, 0.769, 0.082, 0.393, 0.55, 0.255, 0.369, 0.552, 0.452, 0.638, 0.031, 0.592, 0.87, 0.526, 0.747, 0.18, 0.944, 0.975, 0.042, 0.702, 0.849, 0.033, 0.093, 0.594, 0.584, 0.37, 0.848, 0.085, 0.554, 0.301, 0.895, 0.215, 0.259, 0.117, 0.747, 0.828, 0.841, 0.368, 0.985, 0.122, 0.324, 0.696, 0.498, 0.132, 0.878, 0.768, 0.71, 0.812, 0.997, 0.413, 0.277, 0.14, 0.005, 0.596, 0.031, 0.21, 0.405, 0.97, 0.227, 0.885, 0.753, 0.565, 0.555, 0.135, 0.482, 0.536, 0.112, 0.293, 0.516, 0.024, 0.52, 0.78, 0.698, 0.245, 0.772, 0.042, 0.589, 0.631, 0.681, 0.78, 0.993, 0.724, 0.15, 0.392, 0.521, 0.96, 0.596, 0.074, 0.549, 0.193, 0.574, 0.785, 0.818, 0.163, 0.815, 0.169, 0.001, 0.492, 0.381, 0.327, 0.453, 0.737, 0.712, 0.943, 0.426, 0.682, 0.094, 0.79, 0.532, 0.5, 0.9, 0.629, 0.507, 0.142, 0.666, 0.439, 0.391, 0.188, 0.119, 0.782, 0.693, 0.075, 0.276, 0.335, 0.078, 0.969, 0.009, 0.079, 0.422, 0.415, 0.653, 0.058, 0.234, 0.66, 0.565, 0.097, 0.513, 0.345, 0.832, 0.629, 0.43, 0.449, 0.631, 0.521, 0.291, 0.513, 0.719, 0.338, 0.422, 0.844, 0.471, 0.71, 0.363]
global q = [0.843, 0.994, 0.892, 0.887, 0.499, 0.64, 0.222, 0.699, 0.697, 0.974, 0.463, 0.979, 0.967, 0.888, 0.556, 0.921, 0.756, 0.189, 0.917, 0.5, 0.687, 0.72, 0.485, 0.959, 0.984, 0.763, 0.499, 0.435, 0.737, 0.214, 0.591, 0.69, 0.903, 0.884, 0.446, 0.999, 0.1, 0.611, 0.926, 0.84, 0.531, 0.863, 0.565, 0.597, 0.823, 0.7, 0.656, 0.89, 0.963, 0.427, 0.575, 0.871, 0.958, 0.25, 0.92, 0.527, 0.556, 0.766, 0.784, 0.623, 0.417, 0.528, 0.219, 0.338, 0.929, 0.335, 0.933, 0.743, 0.956, 0.866, 0.395, 0.781, 0.461, 0.887, 0.531, 0.666, 0.554, 0.735, 0.456, 0.124, 0.932, 0.866, 0.9, 0.978, 0.703, 0.539, 0.924, 0.989, 0.585, 0.585, 0.737, 0.822, 0.693, 0.932, 0.497, 0.591, 0.945, 0.533, 0.788, 0.676, 0.849, 0.916, 0.396, 0.857, 0.87, 0.847, 0.738, 0.645, 0.69, 0.871, 0.618, 0.988, 0.919, 0.837, 0.894, 0.983, 0.991, 0.937, 0.722, 0.938, 0.96, 0.658, 0.677, 0.592, 0.839, 0.884, 0.597, 0.705, 0.913, 0.913, 0.665, 0.841, 0.421, 0.787, 0.882, 0.869, 0.789, 0.995, 0.849, 0.356, 0.898, 0.779, 0.733, 0.879, 0.882, 0.797, 0.901, 0.998, 0.635, 0.348, 0.864, 0.569, 0.843, 0.902, 0.935, 0.968, 0.989, 0.299, 0.901, 0.916, 0.927, 0.874, 0.951, 0.988, 0.971, 0.403, 0.71, 0.707, 0.853, 0.867, 0.85, 0.961, 0.524, 0.798, 0.946, 0.909, 0.798, 0.962, 0.876, 0.999, 0.962, 0.47, 0.496, 0.633, 0.97, 0.659, 0.739, 0.772, 0.885, 0.649, 0.884, 0.946, 0.962, 0.843, 0.212, 0.949, 0.809, 0.829, 0.742, 0.568, 0.773, 0.912, 0.994, 0.948, 0.733, 0.571, 0.883, 0.809, 0.524, 0.949, 0.875, 0.782, 0.934, 0.947, 0.625, 0.439, 0.66, 0.913, 0.786, 0.986, 0.103, 0.907, 0.767, 0.561, 0.972, 0.303, 0.74, 0.526, 0.475, 0.671, 0.252, 0.9, 0.666, 0.717, 0.369, 0.55, 0.405, 0.882, 0.718, 0.723, 0.46, 0.915, 0.711, 0.598, 0.987, 0.895, 0.389, 0.436, 0.892, 0.976, 0.953, 0.443]
global origin = 1
global destination = 50