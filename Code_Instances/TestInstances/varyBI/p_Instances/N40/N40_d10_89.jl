global arcs = [1 15; 1 25; 1 39; 2 6; 2 20; 2 23; 2 25; 2 32; 3 4; 3 16; 3 20; 3 25; 4 9; 4 23; 4 25; 5 14; 5 20; 5 32; 6 11; 6 16; 6 25; 6 31; 7 16; 7 19; 7 24; 7 28; 7 32; 7 34; 7 35; 7 36; 7 38; 8 3; 8 10; 8 17; 9 7; 9 22; 10 22; 10 27; 10 29; 11 8; 11 18; 11 21; 11 23; 11 26; 11 27; 11 32; 11 40; 12 38; 13 9; 13 28; 13 36; 13 38; 14 2; 14 7; 14 25; 14 27; 15 11; 15 13; 15 18; 15 21; 15 39; 16 11; 16 23; 17 34; 18 13; 18 15; 18 27; 18 32; 19 2; 19 12; 19 13; 19 17; 19 30; 19 31; 19 33; 20 4; 20 6; 20 40; 21 4; 21 10; 21 11; 21 18; 21 20; 21 30; 21 40; 22 2; 22 8; 22 31; 22 34; 22 38; 23 10; 24 5; 24 19; 24 25; 24 28; 24 30; 24 32; 25 3; 25 13; 26 11; 26 14; 26 16; 26 32; 26 34; 27 7; 27 20; 28 16; 28 34; 28 37; 29 8; 29 9; 29 13; 29 18; 29 34; 29 40; 30 5; 30 19; 30 33; 30 34; 31 4; 31 12; 31 23; 31 25; 31 29; 31 33; 32 9; 33 10; 33 15; 33 16; 33 31; 33 38; 34 3; 34 17; 34 19; 34 21; 35 7; 35 24; 35 30; 35 31; 35 34; 35 39; 36 12; 36 13; 36 15; 37 26; 37 38; 38 26; 38 31; 39 15; 39 26; 39 33; 39 36]
global d_x = [1.0, 9.0, 3.0, 3.0, 6.0, 5.0, 3.0, 7.0, 1.0, 9.0, 8.0, 3.0, 10.0, 3.0, 7.0, 5.0, 10.0, 4.0, 2.0, 6.0, 9.0, 9.0, 7.0, 7.0, 6.0, 5.0, 5.0, 4.0, 1.0, 6.0, 10.0, 10.0, 3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 2.0, 7.0, 5.0, 3.0, 6.0, 3.0, 4.0, 6.0, 6.0, 4.0, 10.0, 4.0, 8.0, 6.0, 4.0, 9.0, 8.0, 2.0, 2.0, 9.0, 5.0, 10.0, 9.0, 10.0, 1.0, 8.0, 6.0, 8.0, 1.0, 6.0, 9.0, 1.0, 8.0, 2.0, 9.0, 6.0, 10.0, 10.0, 8.0, 6.0, 9.0, 2.0, 4.0, 5.0, 8.0, 8.0, 5.0, 3.0, 3.0, 7.0, 1.0, 3.0, 10.0, 9.0, 5.0, 5.0, 4.0, 9.0, 2.0, 1.0, 6.0, 4.0, 4.0, 7.0, 7.0, 6.0, 2.0, 3.0, 8.0, 5.0, 2.0, 4.0, 9.0, 1.0, 7.0, 1.0, 10.0, 4.0, 10.0, 6.0, 6.0, 4.0, 8.0, 1.0, 7.0, 5.0, 6.0, 1.0, 4.0, 5.0, 2.0, 5.0, 2.0, 10.0, 7.0, 4.0, 4.0, 4.0, 4.0, 3.0, 10.0, 2.0, 4.0, 4.0, 7.0, 7.0, 5.0, 6.0, 2.0, 4.0, 8.0, 4.0, 2.0, 3.0]
global b_x = 5
global d_y = [5.0, 6.0, 10.0, 10.0, 7.0, 2.0, 2.0, 4.0, 2.0, 6.0, 7.0, 5.0, 5.0, 3.0, 1.0, 5.0, 4.0, 4.0, 8.0, 4.0, 8.0, 2.0, 6.0, 1.0, 5.0, 5.0, 7.0, 9.0, 7.0, 1.0, 2.0, 3.0, 2.0, 8.0, 10.0, 3.0, 10.0, 9.0, 1.0, 1.0, 9.0, 10.0, 7.0, 9.0, 9.0, 9.0, 3.0, 3.0, 3.0, 2.0, 8.0, 10.0, 2.0, 2.0, 4.0, 6.0, 9.0, 6.0, 6.0, 3.0, 10.0, 7.0, 3.0, 2.0, 8.0, 2.0, 5.0, 10.0, 6.0, 10.0, 7.0, 5.0, 5.0, 3.0, 3.0, 5.0, 2.0, 1.0, 9.0, 10.0, 3.0, 7.0, 8.0, 9.0, 6.0, 5.0, 4.0, 3.0, 8.0, 8.0, 8.0, 4.0, 7.0, 6.0, 1.0, 1.0, 5.0, 10.0, 3.0, 9.0, 7.0, 1.0, 9.0, 2.0, 8.0, 9.0, 2.0, 8.0, 3.0, 10.0, 4.0, 8.0, 9.0, 4.0, 5.0, 5.0, 3.0, 5.0, 1.0, 2.0, 1.0, 6.0, 1.0, 8.0, 2.0, 7.0, 10.0, 4.0, 4.0, 1.0, 1.0, 6.0, 3.0, 5.0, 9.0, 8.0, 10.0, 7.0, 1.0, 8.0, 7.0, 1.0, 9.0, 7.0, 10.0, 1.0, 7.0, 8.0, 7.0, 2.0, 4.0, 10.0]
global b_y = 10
global p = [0.937, 0.953, 0.543, 0.82, 0.502, 0.779, 0.715, 0.393, 0.598, 0.666, 0.686, 0.441, 0.968, 0.245, 0.449, 0.58, 0.988, 0.844, 0.29, 0.81, 0.434, 0.871, 0.589, 0.268, 0.05, 0.201, 0.849, 0.448, 0.9, 0.902, 0.881, 0.199, 0.58, 0.527, 0.621, 0.357, 0.733, 0.24, 0.064, 0.7, 0.59, 0.769, 0.091, 0.711, 0.003, 0.687, 0.052, 0.372, 0.283, 0.049, 0.198, 0.903, 0.783, 0.753, 0.041, 0.811, 0.328, 0.116, 0.902, 0.84, 0.775, 0.773, 0.874, 0.505, 0.846, 0.025, 0.499, 0.968, 0.356, 0.507, 0.795, 0.847, 0.545, 0.04, 0.332, 0.149, 0.361, 0.59, 0.007, 0.08, 0.819, 0.74, 0.723, 0.427, 0.024, 0.409, 0.896, 0.126, 0.158, 0.233, 0.581, 0.371, 0.888, 0.53, 0.699, 0.444, 0.676, 0.285, 0.291, 0.142, 0.26, 0.177, 0.938, 0.107, 0.915, 0.082, 0.813, 0.059, 0.116, 0.531, 0.038, 0.403, 0.978, 0.017, 0.489, 0.914, 0.627, 0.439, 0.651, 0.736, 0.143, 0.638, 0.119, 0.687, 0.124, 0.562, 0.8, 0.553, 0.67, 0.764, 0.806, 0.768, 0.345, 0.301, 0.077, 0.73, 0.043, 0.581, 0.425, 0.889, 0.217, 0.636, 0.062, 0.848, 0.267, 0.49, 0.53, 0.014, 0.783, 0.256, 0.236, 0.052]
global q = [0.942, 0.966, 0.804, 0.865, 0.984, 0.954, 0.831, 0.885, 0.944, 0.971, 0.971, 0.864, 0.97, 0.994, 0.822, 0.863, 0.991, 0.903, 0.663, 0.845, 0.463, 0.887, 0.934, 0.93, 0.798, 0.642, 0.972, 0.66, 0.973, 0.927, 0.947, 0.533, 0.673, 0.853, 0.978, 0.967, 0.791, 0.445, 0.867, 0.707, 0.826, 0.95, 0.362, 0.776, 0.469, 0.714, 0.575, 0.379, 0.615, 0.843, 0.334, 0.994, 0.8, 0.839, 0.887, 0.834, 0.674, 0.572, 0.933, 0.94, 0.903, 0.978, 0.964, 0.61, 0.993, 0.17, 0.608, 0.987, 0.958, 0.64, 0.97, 0.997, 0.557, 0.712, 0.979, 0.83, 0.916, 0.756, 0.922, 0.532, 0.896, 0.844, 0.975, 0.754, 0.749, 0.874, 0.898, 0.778, 0.983, 0.946, 0.794, 0.836, 0.928, 0.793, 0.923, 0.591, 0.855, 0.668, 0.634, 0.277, 0.512, 0.95, 0.989, 0.394, 0.992, 0.437, 0.95, 0.576, 0.395, 0.734, 0.578, 0.836, 0.999, 0.854, 0.862, 0.984, 0.782, 0.501, 0.679, 0.846, 0.813, 0.936, 0.501, 0.98, 0.366, 0.62, 0.852, 0.624, 0.871, 0.77, 0.975, 0.799, 0.46, 0.908, 0.253, 0.799, 0.673, 0.662, 0.646, 0.902, 0.534, 0.881, 0.089, 0.959, 0.553, 0.564, 0.645, 0.667, 0.925, 0.492, 0.875, 0.316]
global origin = 1
global destination = 40
