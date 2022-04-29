global arcs = [1 12; 1 21; 1 22; 1 29; 2 18; 2 26; 3 12; 3 20; 4 8; 4 12; 4 15; 4 16; 4 18; 4 24; 5 10; 5 21; 6 5; 6 16; 7 9; 7 23; 7 24; 8 3; 8 9; 9 4; 9 19; 9 20; 10 17; 10 23; 11 2; 11 3; 11 6; 11 9; 11 27; 11 29; 12 5; 12 25; 13 9; 13 28; 14 10; 14 11; 14 12; 14 19; 15 2; 15 13; 15 14; 15 18; 15 21; 15 24; 15 25; 15 26; 16 19; 16 24; 17 10; 17 11; 17 12; 17 20; 17 24; 18 13; 18 30; 19 10; 19 13; 19 28; 20 11; 20 24; 21 3; 21 4; 21 5; 21 14; 22 11; 23 15; 23 22; 23 28; 23 30; 24 5; 24 6; 24 17; 24 20; 24 29; 25 21; 26 6; 26 7; 26 16; 26 25; 26 30; 27 26; 28 3; 28 8; 28 27; 29 3; 29 4; 29 6]
global d_x = [2.0, 5.0, 2.0, 5.0, 10.0, 3.0, 2.0, 5.0, 6.0, 5.0, 10.0, 2.0, 3.0, 8.0, 3.0, 9.0, 10.0, 7.0, 10.0, 5.0, 10.0, 7.0, 1.0, 2.0, 7.0, 10.0, 10.0, 2.0, 4.0, 5.0, 8.0, 3.0, 6.0, 5.0, 4.0, 1.0, 7.0, 1.0, 2.0, 8.0, 1.0, 5.0, 2.0, 4.0, 2.0, 6.0, 1.0, 9.0, 2.0, 1.0, 7.0, 3.0, 2.0, 2.0, 6.0, 6.0, 1.0, 1.0, 6.0, 3.0, 10.0, 4.0, 4.0, 6.0, 10.0, 2.0, 10.0, 10.0, 5.0, 2.0, 7.0, 8.0, 6.0, 7.0, 6.0, 10.0, 2.0, 10.0, 3.0, 4.0, 10.0, 4.0, 1.0, 7.0, 8.0, 9.0, 6.0, 3.0, 5.0, 10.0, 6.0]
global b_x = 5
global d_y = [10.0, 5.0, 3.0, 1.0, 7.0, 9.0, 2.0, 7.0, 2.0, 9.0, 2.0, 5.0, 10.0, 4.0, 7.0, 4.0, 10.0, 2.0, 9.0, 1.0, 1.0, 2.0, 2.0, 1.0, 4.0, 7.0, 3.0, 3.0, 5.0, 6.0, 6.0, 6.0, 10.0, 7.0, 6.0, 2.0, 4.0, 2.0, 5.0, 8.0, 10.0, 9.0, 3.0, 2.0, 3.0, 7.0, 2.0, 8.0, 7.0, 4.0, 7.0, 8.0, 7.0, 1.0, 3.0, 3.0, 7.0, 4.0, 10.0, 1.0, 1.0, 1.0, 10.0, 4.0, 4.0, 8.0, 8.0, 7.0, 6.0, 8.0, 3.0, 3.0, 3.0, 6.0, 2.0, 4.0, 7.0, 4.0, 1.0, 4.0, 3.0, 7.0, 10.0, 2.0, 3.0, 3.0, 9.0, 3.0, 9.0, 4.0, 5.0]
global b_y = 10
global p = [0.593, 0.91, 0.367, 0.813, 0.623, 0.207, 0.392, 0.875, 0.838, 0.37, 0.019, 0.327, 0.895, 0.987, 0.857, 0.792, 0.677, 0.511, 0.751, 0.845, 0.661, 0.406, 0.51, 0.232, 0.574, 0.427, 0.956, 0.946, 0.002, 0.115, 0.551, 0.89, 0.362, 0.745, 0.753, 0.882, 0.215, 0.67, 0.771, 0.601, 0.661, 0.819, 0.711, 0.633, 0.359, 0.954, 0.697, 0.176, 0.141, 0.642, 0.437, 0.168, 0.508, 0.568, 0.367, 0.773, 0.626, 0.236, 0.531, 0.489, 0.836, 0.748, 0.087, 0.852, 0.545, 0.834, 0.281, 0.532, 0.583, 0.578, 0.777, 0.497, 0.428, 0.987, 0.145, 0.33, 0.976, 0.29, 0.532, 0.062, 0.226, 0.811, 0.294, 0.922, 0.175, 0.761, 0.579, 0.743, 0.335, 0.501, 0.918]
global q = [0.642, 0.949, 0.438, 0.946, 0.87, 0.394, 0.762, 0.966, 0.874, 0.989, 0.939, 0.733, 0.923, 0.993, 0.862, 0.862, 0.869, 0.837, 0.945, 0.856, 0.787, 0.859, 0.709, 0.523, 0.828, 0.617, 0.99, 0.946, 0.073, 0.185, 0.994, 0.917, 0.494, 0.871, 0.788, 0.936, 0.516, 0.713, 0.837, 0.829, 0.716, 0.928, 0.837, 0.69, 0.928, 0.963, 0.992, 0.79, 0.38, 0.716, 0.93, 0.216, 0.662, 0.896, 0.951, 0.881, 0.837, 0.825, 0.884, 0.646, 0.889, 0.827, 0.823, 0.896, 0.852, 0.868, 0.623, 0.957, 0.891, 0.899, 0.932, 0.561, 0.88, 0.998, 0.996, 0.985, 0.991, 0.298, 0.592, 0.577, 0.391, 0.905, 0.821, 0.938, 0.374, 0.782, 0.761, 0.869, 0.521, 0.72, 0.997]
global origin = 1
global destination = 30
