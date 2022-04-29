global arcs = [1 31; 2 11; 2 12; 2 16; 2 17; 2 22; 2 24; 3 15; 3 24; 3 33; 3 34; 4 8; 4 10; 4 11; 4 18; 4 21; 5 3; 5 14; 5 15; 6 4; 6 28; 6 30; 6 34; 7 19; 7 24; 7 28; 8 26; 8 35; 9 8; 9 16; 9 17; 9 23; 10 17; 10 28; 11 2; 11 3; 12 28; 13 5; 13 10; 14 6; 14 8; 14 9; 14 13; 14 30; 15 8; 15 9; 15 22; 16 3; 16 26; 16 35; 17 7; 17 15; 18 8; 18 14; 18 27; 19 9; 19 12; 19 22; 19 24; 19 25; 19 32; 20 8; 20 10; 20 13; 20 14; 20 18; 20 23; 20 24; 20 30; 20 31; 21 2; 21 5; 21 7; 21 30; 22 2; 22 7; 22 14; 22 24; 22 30; 22 33; 23 2; 23 6; 23 11; 23 13; 23 30; 23 33; 24 7; 24 21; 24 33; 25 4; 25 8; 25 29; 25 33; 26 10; 26 17; 26 22; 26 25; 26 35; 27 5; 27 6; 27 9; 28 12; 28 14; 28 21; 28 22; 28 30; 28 32; 28 34; 29 8; 30 4; 30 17; 30 19; 30 20; 31 6; 31 8; 31 17; 31 26; 31 30; 32 16; 32 25; 33 3; 33 6; 33 14; 33 29; 33 30; 34 28]
global d_x = [7.0, 9.0, 6.0, 7.0, 10.0, 6.0, 7.0, 9.0, 2.0, 5.0, 10.0, 8.0, 4.0, 7.0, 10.0, 6.0, 9.0, 10.0, 2.0, 10.0, 8.0, 9.0, 10.0, 6.0, 4.0, 9.0, 8.0, 9.0, 6.0, 8.0, 2.0, 3.0, 9.0, 4.0, 3.0, 5.0, 10.0, 6.0, 6.0, 1.0, 5.0, 10.0, 10.0, 8.0, 8.0, 10.0, 1.0, 4.0, 10.0, 5.0, 6.0, 10.0, 1.0, 9.0, 6.0, 9.0, 8.0, 5.0, 8.0, 6.0, 3.0, 10.0, 3.0, 3.0, 8.0, 10.0, 6.0, 7.0, 6.0, 3.0, 5.0, 5.0, 5.0, 3.0, 8.0, 3.0, 2.0, 2.0, 6.0, 9.0, 8.0, 9.0, 9.0, 8.0, 1.0, 5.0, 2.0, 6.0, 5.0, 6.0, 4.0, 1.0, 7.0, 3.0, 6.0, 3.0, 6.0, 5.0, 7.0, 1.0, 6.0, 9.0, 10.0, 8.0, 4.0, 5.0, 5.0, 8.0, 2.0, 10.0, 2.0, 5.0, 7.0, 5.0, 7.0, 3.0, 8.0, 9.0, 4.0, 6.0, 7.0, 10.0, 4.0, 6.0, 7.0, 8.0]
global b_x = 5
global d_y = [1.0, 5.0, 7.0, 6.0, 4.0, 5.0, 10.0, 7.0, 4.0, 7.0, 3.0, 1.0, 3.0, 4.0, 2.0, 4.0, 7.0, 1.0, 1.0, 7.0, 6.0, 9.0, 10.0, 8.0, 8.0, 8.0, 3.0, 6.0, 7.0, 7.0, 10.0, 7.0, 1.0, 4.0, 9.0, 1.0, 7.0, 10.0, 4.0, 1.0, 1.0, 4.0, 7.0, 10.0, 9.0, 10.0, 9.0, 4.0, 2.0, 4.0, 7.0, 8.0, 3.0, 2.0, 2.0, 4.0, 1.0, 1.0, 6.0, 8.0, 2.0, 4.0, 10.0, 3.0, 5.0, 8.0, 9.0, 10.0, 2.0, 2.0, 3.0, 5.0, 3.0, 7.0, 10.0, 10.0, 2.0, 7.0, 6.0, 10.0, 3.0, 7.0, 10.0, 4.0, 7.0, 3.0, 3.0, 2.0, 9.0, 10.0, 9.0, 8.0, 4.0, 7.0, 10.0, 1.0, 7.0, 9.0, 5.0, 2.0, 3.0, 3.0, 10.0, 2.0, 6.0, 3.0, 1.0, 1.0, 1.0, 1.0, 6.0, 8.0, 4.0, 7.0, 5.0, 3.0, 6.0, 1.0, 7.0, 6.0, 5.0, 5.0, 4.0, 10.0, 6.0, 2.0]
global b_y = 10
global p = [0.918, 0.255, 0.546, 0.783, 0.71, 0.885, 0.577, 0.654, 0.164, 0.57, 0.834, 0.115, 0.763, 0.307, 0.505, 0.852, 0.346, 0.588, 0.156, 0.803, 0.813, 0.228, 0.172, 0.239, 0.295, 0.327, 0.583, 0.491, 0.861, 0.149, 0.832, 0.2, 0.395, 0.436, 0.726, 0.947, 0.309, 0.164, 0.119, 0.206, 0.682, 0.557, 0.458, 0.266, 0.691, 0.926, 0.463, 0.685, 0.676, 0.793, 0.328, 0.124, 0.192, 0.073, 0.295, 0.383, 0.985, 0.108, 0.853, 0.844, 0.614, 0.498, 0.853, 0.509, 0.929, 0.165, 0.152, 0.73, 0.58, 0.354, 0.034, 0.689, 0.008, 0.06, 0.04, 0.323, 0.727, 0.609, 0.987, 0.774, 0.37, 0.112, 0.829, 0.261, 0.397, 0.417, 0.925, 0.193, 0.374, 0.611, 0.812, 0.925, 0.535, 0.654, 0.03, 0.663, 0.006, 0.428, 0.937, 0.896, 0.119, 0.564, 0.378, 0.666, 0.889, 0.331, 0.7, 0.101, 0.079, 0.675, 0.862, 0.11, 0.348, 0.879, 0.686, 0.436, 0.867, 0.649, 0.227, 0.565, 0.994, 0.51, 0.547, 0.815, 0.718, 0.372]
global q = [0.939, 0.299, 0.778, 0.81, 0.789, 0.956, 0.916, 0.711, 0.857, 0.879, 0.982, 0.335, 0.825, 0.949, 0.96, 0.961, 0.567, 0.859, 0.717, 0.99, 0.815, 0.437, 0.524, 0.358, 0.304, 0.84, 0.725, 0.944, 0.919, 0.498, 0.973, 0.462, 0.639, 0.901, 0.807, 0.962, 0.67, 0.307, 0.37, 0.814, 0.918, 0.596, 0.475, 0.703, 0.886, 0.947, 0.812, 0.696, 0.915, 0.991, 0.778, 0.62, 0.515, 0.406, 0.492, 0.575, 0.987, 0.813, 0.895, 0.877, 0.708, 0.804, 0.908, 0.628, 0.939, 0.914, 0.903, 0.875, 0.624, 0.756, 0.283, 0.897, 0.966, 0.722, 0.675, 0.329, 0.935, 0.834, 0.991, 0.866, 0.842, 0.223, 0.847, 0.938, 0.664, 0.964, 0.967, 0.27, 0.873, 0.67, 0.94, 0.952, 0.968, 0.931, 0.116, 0.757, 0.837, 0.487, 0.974, 0.949, 0.469, 0.63, 0.607, 0.806, 0.949, 0.56, 0.886, 0.134, 0.958, 0.864, 0.934, 0.439, 0.414, 0.978, 0.888, 0.688, 0.964, 0.695, 0.502, 0.827, 0.997, 0.795, 0.603, 0.951, 0.742, 0.586]
global origin = 1
global destination = 35