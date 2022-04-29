global arcs = [1 2; 1 13; 1 24; 1 25; 1 40; 2 30; 2 31; 3 16; 3 31; 4 2; 4 23; 5 2; 5 3; 5 8; 5 36; 6 13; 6 19; 6 21; 6 27; 6 35; 7 2; 7 6; 7 9; 7 20; 7 33; 7 38; 8 14; 8 21; 8 28; 8 31; 9 26; 9 30; 9 36; 10 5; 10 8; 11 6; 11 7; 11 30; 11 40; 12 3; 12 8; 12 10; 12 13; 12 26; 12 30; 12 35; 12 37; 12 40; 13 6; 13 31; 13 34; 14 26; 14 32; 15 29; 15 37; 15 39; 15 40; 16 13; 16 15; 16 38; 17 16; 17 19; 17 24; 17 34; 18 16; 18 29; 18 30; 19 29; 19 34; 19 37; 20 2; 20 17; 20 32; 20 33; 21 6; 21 8; 21 10; 21 13; 21 19; 21 34; 22 11; 22 18; 22 21; 22 35; 23 2; 23 6; 23 17; 23 19; 23 29; 24 10; 24 15; 24 25; 25 3; 25 13; 25 20; 25 21; 25 33; 25 39; 26 33; 27 16; 27 17; 27 22; 27 26; 28 4; 28 19; 28 26; 28 36; 29 13; 29 26; 30 4; 30 11; 30 12; 30 23; 30 24; 30 25; 30 38; 31 9; 31 12; 31 14; 31 30; 31 36; 31 40; 32 18; 32 24; 32 29; 33 22; 33 32; 33 36; 33 37; 34 3; 34 10; 34 13; 34 20; 34 23; 34 39; 35 21; 35 24; 35 40; 36 7; 36 8; 36 16; 36 20; 36 23; 37 4; 37 16; 37 18; 37 23; 38 5; 38 11; 38 26; 39 34]
global d_x = [2.0, 5.0, 5.0, 9.0, 6.0, 4.0, 4.0, 5.0, 9.0, 1.0, 4.0, 4.0, 9.0, 7.0, 10.0, 2.0, 4.0, 4.0, 3.0, 4.0, 8.0, 6.0, 4.0, 1.0, 9.0, 5.0, 5.0, 2.0, 7.0, 8.0, 7.0, 10.0, 2.0, 6.0, 4.0, 7.0, 3.0, 8.0, 4.0, 8.0, 1.0, 2.0, 4.0, 10.0, 1.0, 8.0, 3.0, 7.0, 3.0, 1.0, 6.0, 2.0, 9.0, 6.0, 3.0, 4.0, 4.0, 6.0, 9.0, 8.0, 8.0, 5.0, 1.0, 2.0, 6.0, 10.0, 5.0, 5.0, 4.0, 7.0, 1.0, 2.0, 2.0, 10.0, 3.0, 3.0, 4.0, 7.0, 5.0, 10.0, 10.0, 10.0, 3.0, 7.0, 9.0, 10.0, 3.0, 7.0, 9.0, 10.0, 6.0, 1.0, 6.0, 9.0, 1.0, 10.0, 4.0, 8.0, 5.0, 4.0, 9.0, 9.0, 9.0, 2.0, 7.0, 4.0, 4.0, 5.0, 6.0, 2.0, 10.0, 1.0, 1.0, 7.0, 7.0, 9.0, 7.0, 7.0, 7.0, 5.0, 1.0, 10.0, 3.0, 2.0, 4.0, 6.0, 2.0, 1.0, 7.0, 1.0, 10.0, 8.0, 8.0, 9.0, 1.0, 5.0, 8.0, 7.0, 6.0, 4.0, 8.0, 6.0, 5.0, 2.0, 9.0, 7.0, 10.0, 9.0, 3.0, 9.0, 10.0]
global b_x = 5
global d_y = [4.0, 3.0, 2.0, 1.0, 2.0, 2.0, 7.0, 5.0, 6.0, 8.0, 3.0, 7.0, 4.0, 10.0, 9.0, 8.0, 5.0, 3.0, 6.0, 7.0, 10.0, 7.0, 1.0, 4.0, 9.0, 10.0, 8.0, 3.0, 6.0, 9.0, 9.0, 8.0, 1.0, 10.0, 1.0, 2.0, 4.0, 9.0, 3.0, 8.0, 5.0, 5.0, 8.0, 3.0, 7.0, 8.0, 3.0, 4.0, 6.0, 4.0, 6.0, 5.0, 10.0, 5.0, 10.0, 3.0, 7.0, 7.0, 10.0, 8.0, 8.0, 3.0, 9.0, 8.0, 10.0, 2.0, 1.0, 9.0, 9.0, 9.0, 9.0, 7.0, 6.0, 9.0, 4.0, 5.0, 4.0, 2.0, 1.0, 3.0, 1.0, 3.0, 6.0, 9.0, 5.0, 10.0, 3.0, 8.0, 6.0, 3.0, 1.0, 9.0, 8.0, 8.0, 1.0, 6.0, 2.0, 4.0, 9.0, 6.0, 3.0, 10.0, 7.0, 6.0, 9.0, 9.0, 4.0, 5.0, 6.0, 2.0, 3.0, 1.0, 10.0, 2.0, 10.0, 2.0, 3.0, 7.0, 4.0, 1.0, 9.0, 5.0, 8.0, 1.0, 10.0, 8.0, 6.0, 2.0, 7.0, 4.0, 7.0, 1.0, 2.0, 2.0, 9.0, 5.0, 9.0, 3.0, 6.0, 8.0, 4.0, 6.0, 9.0, 7.0, 4.0, 5.0, 9.0, 5.0, 2.0, 5.0, 1.0]
global b_y = 10
global p = [0.836, 0.978, 0.702, 0.14, 0.919, 0.359, 0.81, 0.156, 0.492, 0.597, 0.743, 0.878, 0.74, 0.257, 0.436, 0.963, 0.02, 0.455, 0.437, 0.668, 0.49, 0.233, 0.688, 0.486, 0.961, 0.014, 0.727, 0.339, 0.137, 0.027, 0.066, 0.334, 0.591, 0.478, 0.573, 0.213, 0.904, 0.961, 0.051, 0.378, 0.932, 0.041, 0.183, 0.587, 0.057, 0.614, 0.361, 0.386, 0.277, 0.613, 0.273, 0.159, 0.074, 0.878, 0.604, 0.828, 0.934, 0.36, 0.472, 0.882, 0.078, 0.606, 0.22, 0.159, 0.138, 0.69, 0.716, 0.396, 0.671, 0.666, 0.565, 0.152, 0.142, 0.93, 0.148, 0.744, 0.169, 0.42, 0.696, 0.642, 0.695, 0.978, 0.561, 0.324, 0.581, 0.541, 0.194, 0.372, 0.364, 0.403, 0.217, 0.485, 0.856, 0.43, 0.087, 0.34, 0.041, 0.426, 0.244, 0.871, 0.94, 0.362, 0.803, 0.222, 0.038, 0.009, 0.538, 0.99, 0.705, 0.177, 0.824, 0.804, 0.329, 0.601, 0.28, 0.206, 0.53, 0.041, 0.234, 0.684, 0.287, 0.706, 0.214, 0.596, 0.911, 0.198, 0.721, 0.252, 0.616, 0.335, 0.128, 0.894, 0.149, 0.14, 0.268, 0.249, 0.256, 0.435, 0.571, 0.71, 0.145, 0.745, 0.688, 0.246, 0.69, 0.664, 0.883, 0.547, 0.533, 0.74, 0.518]
global q = [0.917, 0.986, 0.83, 0.355, 0.949, 0.617, 0.872, 0.23, 0.527, 0.798, 0.875, 0.922, 0.774, 0.295, 0.926, 0.995, 0.364, 0.883, 0.582, 0.777, 0.856, 0.396, 0.8, 0.972, 0.968, 0.757, 0.854, 0.679, 0.29, 0.848, 0.895, 0.837, 0.834, 0.919, 0.632, 0.561, 0.972, 0.968, 0.149, 0.689, 0.994, 0.478, 0.557, 0.824, 0.161, 0.896, 0.795, 0.97, 0.378, 0.966, 0.624, 0.579, 0.26, 0.983, 0.753, 0.85, 0.974, 0.59, 0.867, 0.91, 0.41, 0.795, 0.338, 0.534, 0.173, 0.77, 0.964, 0.668, 0.83, 0.724, 0.631, 0.262, 0.917, 0.994, 0.652, 0.789, 0.875, 0.918, 0.992, 0.968, 0.807, 0.994, 0.953, 0.515, 0.645, 0.573, 0.209, 0.414, 0.52, 0.782, 0.567, 0.822, 0.857, 0.449, 0.202, 0.595, 0.586, 0.546, 0.919, 0.946, 0.96, 0.79, 0.812, 0.418, 0.457, 0.296, 0.894, 0.993, 0.753, 0.229, 0.866, 0.894, 0.527, 0.663, 0.382, 0.403, 0.615, 0.089, 0.435, 0.721, 0.454, 0.964, 0.454, 0.767, 0.973, 0.27, 0.915, 0.798, 0.841, 0.899, 0.205, 0.966, 0.318, 0.347, 0.364, 0.614, 0.472, 0.546, 0.916, 0.856, 0.365, 0.942, 0.953, 0.651, 0.862, 0.829, 0.976, 0.726, 0.666, 0.79, 0.843]
global origin = 1
global destination = 40