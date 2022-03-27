global arcs = [1 3; 1 27; 1 32; 2 10; 2 16; 2 17; 2 24; 2 40; 3 4; 3 9; 3 18; 3 22; 3 28; 3 29; 3 40; 4 7; 4 26; 5 3; 5 16; 5 26; 6 5; 6 10; 6 23; 7 4; 7 10; 7 12; 7 20; 7 28; 7 34; 8 10; 8 36; 9 6; 9 10; 9 19; 10 2; 10 39; 11 10; 11 27; 12 4; 12 10; 12 21; 12 38; 13 3; 13 9; 14 6; 15 5; 15 21; 16 28; 17 18; 17 35; 18 15; 18 22; 18 23; 18 27; 19 5; 19 6; 19 7; 19 25; 19 34; 19 35; 20 12; 20 16; 20 37; 21 6; 21 7; 21 13; 21 30; 21 31; 21 39; 22 2; 22 3; 22 6; 23 14; 23 26; 23 31; 23 39; 24 6; 24 7; 24 18; 24 19; 24 29; 25 11; 25 35; 25 36; 26 5; 26 13; 26 18; 26 24; 26 33; 26 38; 27 5; 27 8; 27 24; 27 29; 28 25; 29 37; 29 40; 30 10; 30 16; 31 16; 31 24; 31 25; 31 26; 31 39; 32 2; 32 4; 32 16; 32 21; 33 7; 33 15; 33 32; 34 4; 34 14; 34 22; 34 24; 34 29; 35 2; 35 11; 35 14; 35 16; 35 18; 35 23; 35 27; 35 30; 35 32; 35 34; 35 36; 36 4; 36 8; 36 15; 36 24; 36 31; 37 2; 37 4; 37 7; 37 9; 37 29; 37 34; 37 39; 38 4; 38 29; 38 39; 39 4; 39 7; 39 28]
global d_x = [2.0, 10.0, 7.0, 5.0, 8.0, 8.0, 8.0, 7.0, 4.0, 5.0, 7.0, 10.0, 2.0, 5.0, 5.0, 8.0, 6.0, 6.0, 6.0, 1.0, 1.0, 4.0, 2.0, 3.0, 4.0, 4.0, 2.0, 9.0, 3.0, 3.0, 4.0, 10.0, 9.0, 5.0, 4.0, 6.0, 1.0, 2.0, 1.0, 5.0, 10.0, 5.0, 1.0, 6.0, 10.0, 8.0, 3.0, 10.0, 2.0, 2.0, 2.0, 4.0, 10.0, 6.0, 9.0, 4.0, 1.0, 10.0, 6.0, 2.0, 6.0, 6.0, 4.0, 10.0, 10.0, 3.0, 1.0, 8.0, 5.0, 9.0, 5.0, 10.0, 1.0, 3.0, 4.0, 3.0, 6.0, 4.0, 1.0, 5.0, 9.0, 5.0, 6.0, 5.0, 9.0, 5.0, 2.0, 8.0, 8.0, 5.0, 5.0, 5.0, 9.0, 6.0, 8.0, 1.0, 2.0, 6.0, 6.0, 3.0, 6.0, 2.0, 1.0, 9.0, 1.0, 9.0, 5.0, 8.0, 3.0, 9.0, 10.0, 8.0, 10.0, 2.0, 1.0, 6.0, 5.0, 9.0, 1.0, 9.0, 7.0, 3.0, 9.0, 3.0, 2.0, 2.0, 7.0, 1.0, 1.0, 4.0, 6.0, 10.0, 6.0, 1.0, 6.0, 3.0, 6.0, 3.0, 10.0, 4.0, 10.0, 4.0, 7.0, 10.0, 8.0]
global b_x = 5
global d_y = [2.0, 5.0, 5.0, 7.0, 4.0, 7.0, 3.0, 3.0, 2.0, 3.0, 6.0, 7.0, 3.0, 6.0, 5.0, 10.0, 5.0, 1.0, 4.0, 7.0, 10.0, 10.0, 8.0, 10.0, 1.0, 8.0, 5.0, 3.0, 1.0, 7.0, 10.0, 3.0, 3.0, 7.0, 10.0, 7.0, 6.0, 5.0, 5.0, 8.0, 10.0, 2.0, 9.0, 4.0, 8.0, 1.0, 5.0, 1.0, 2.0, 8.0, 7.0, 3.0, 5.0, 4.0, 10.0, 2.0, 10.0, 7.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 3.0, 3.0, 6.0, 8.0, 9.0, 7.0, 1.0, 9.0, 4.0, 4.0, 6.0, 4.0, 4.0, 10.0, 1.0, 3.0, 3.0, 6.0, 9.0, 3.0, 6.0, 1.0, 3.0, 4.0, 6.0, 2.0, 10.0, 7.0, 7.0, 6.0, 6.0, 6.0, 5.0, 9.0, 5.0, 7.0, 10.0, 8.0, 5.0, 9.0, 7.0, 8.0, 4.0, 5.0, 2.0, 1.0, 9.0, 10.0, 1.0, 9.0, 7.0, 9.0, 8.0, 6.0, 9.0, 8.0, 7.0, 8.0, 6.0, 10.0, 1.0, 2.0, 4.0, 9.0, 3.0, 7.0, 8.0, 2.0, 7.0, 10.0, 8.0, 10.0, 3.0, 6.0, 8.0, 3.0, 2.0, 3.0, 10.0, 9.0, 1.0]
global b_y = 10
global p = [0.411, 0.826, 0.096, 0.978, 0.628, 0.611, 0.615, 0.921, 0.516, 0.256, 0.919, 0.319, 0.058, 0.338, 0.625, 0.828, 0.734, 0.355, 0.812, 0.28, 0.237, 0.975, 0.078, 0.107, 0.317, 0.326, 0.16, 0.977, 0.253, 0.01, 0.111, 0.858, 0.066, 0.349, 0.095, 0.406, 0.01, 0.623, 0.254, 0.663, 0.694, 0.549, 0.438, 0.849, 0.251, 0.313, 0.53, 0.314, 0.505, 0.641, 0.354, 0.596, 0.403, 0.485, 0.811, 0.504, 0.968, 0.539, 0.707, 0.854, 0.192, 0.323, 0.113, 0.007, 0.126, 0.531, 0.774, 0.478, 0.44, 0.868, 0.876, 0.846, 0.576, 0.24, 0.094, 0.779, 0.96, 0.259, 0.24, 0.691, 0.598, 0.637, 0.857, 0.356, 0.512, 0.058, 0.856, 0.757, 0.909, 0.792, 0.899, 0.374, 0.353, 0.601, 0.746, 0.959, 0.282, 0.653, 0.548, 0.298, 0.939, 0.522, 0.684, 0.254, 0.427, 0.677, 0.303, 0.609, 0.326, 0.961, 0.81, 0.73, 0.29, 0.521, 0.335, 0.452, 0.19, 0.138, 0.836, 0.325, 0.038, 0.112, 0.529, 0.368, 0.689, 0.439, 0.282, 0.824, 0.736, 0.582, 0.187, 0.737, 0.914, 0.011, 0.105, 0.166, 0.205, 0.036, 0.482, 0.567, 0.037, 0.376, 0.957, 0.72, 0.79]
global q = [0.856, 0.868, 0.936, 0.999, 0.931, 0.669, 0.678, 0.955, 0.958, 0.623, 0.96, 0.483, 0.997, 0.983, 0.661, 0.872, 0.748, 0.508, 0.904, 0.335, 0.435, 0.977, 0.302, 0.543, 0.767, 0.975, 0.343, 0.995, 0.923, 0.579, 0.958, 0.991, 0.076, 0.664, 0.577, 0.484, 0.136, 0.808, 0.722, 0.745, 0.99, 0.695, 0.567, 0.947, 0.731, 0.392, 0.678, 0.44, 0.595, 0.738, 0.447, 0.959, 0.865, 0.651, 0.961, 0.784, 0.991, 0.796, 0.974, 0.93, 0.659, 0.354, 0.617, 0.546, 0.598, 0.572, 0.947, 0.743, 0.454, 0.986, 0.972, 0.974, 0.708, 0.96, 0.474, 0.791, 0.989, 0.622, 0.73, 0.726, 0.633, 0.68, 0.95, 0.778, 0.57, 0.972, 0.917, 0.842, 0.937, 0.882, 0.927, 0.846, 0.47, 0.765, 0.972, 0.981, 0.284, 0.771, 0.777, 0.472, 0.952, 0.773, 0.849, 0.51, 0.974, 0.865, 0.583, 0.802, 0.589, 0.982, 0.812, 0.731, 0.857, 0.622, 0.673, 0.711, 0.44, 0.781, 0.871, 0.643, 0.591, 0.303, 0.619, 0.649, 0.941, 0.872, 0.549, 0.955, 0.844, 0.741, 0.214, 0.847, 0.994, 0.201, 0.372, 0.959, 0.236, 0.189, 0.943, 0.626, 0.051, 0.634, 0.972, 0.838, 0.931]
global origin = 1
global destination = 40