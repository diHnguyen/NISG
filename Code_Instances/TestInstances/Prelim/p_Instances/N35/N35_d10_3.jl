global arcs = [1 18; 1 20; 2 22; 2 25; 3 9; 3 14; 3 16; 3 26; 4 11; 4 12; 4 17; 5 14; 5 23; 5 27; 5 28; 5 33; 5 35; 6 13; 6 21; 6 30; 7 24; 8 24; 9 10; 9 28; 9 35; 10 2; 10 5; 10 6; 10 8; 10 11; 10 14; 10 16; 10 22; 10 23; 10 24; 10 26; 10 27; 10 29; 10 32; 11 2; 11 13; 12 16; 12 24; 12 34; 13 11; 13 23; 13 29; 13 30; 13 33; 14 31; 14 32; 14 34; 15 3; 15 4; 15 10; 15 29; 15 31; 16 2; 16 4; 16 26; 16 27; 16 29; 17 8; 17 18; 17 25; 18 23; 18 28; 19 20; 19 27; 19 31; 20 17; 20 34; 21 8; 21 15; 21 26; 21 27; 21 29; 21 34; 21 35; 22 29; 23 12; 23 13; 23 14; 23 16; 23 26; 23 28; 23 29; 24 2; 24 19; 24 21; 24 26; 24 35; 25 2; 25 9; 25 28; 26 3; 26 6; 26 29; 27 2; 27 10; 27 24; 28 7; 28 24; 29 13; 29 20; 29 28; 29 33; 30 14; 30 25; 31 17; 31 22; 31 33; 32 13; 32 19; 32 22; 32 35; 33 3; 33 17; 33 21; 34 4; 34 32]
global d_x = [10.0, 9.0, 1.0, 8.0, 1.0, 3.0, 8.0, 4.0, 2.0, 2.0, 1.0, 10.0, 1.0, 9.0, 10.0, 8.0, 8.0, 5.0, 3.0, 8.0, 10.0, 9.0, 9.0, 6.0, 9.0, 7.0, 9.0, 2.0, 6.0, 10.0, 5.0, 2.0, 5.0, 7.0, 10.0, 8.0, 8.0, 1.0, 5.0, 9.0, 8.0, 10.0, 4.0, 6.0, 10.0, 4.0, 1.0, 8.0, 6.0, 5.0, 2.0, 2.0, 10.0, 7.0, 4.0, 4.0, 3.0, 1.0, 1.0, 7.0, 4.0, 9.0, 4.0, 6.0, 3.0, 7.0, 6.0, 5.0, 1.0, 2.0, 1.0, 4.0, 3.0, 8.0, 4.0, 3.0, 10.0, 6.0, 3.0, 10.0, 6.0, 10.0, 4.0, 3.0, 4.0, 7.0, 6.0, 6.0, 10.0, 8.0, 3.0, 8.0, 4.0, 3.0, 6.0, 3.0, 8.0, 9.0, 5.0, 10.0, 7.0, 7.0, 4.0, 8.0, 4.0, 10.0, 7.0, 4.0, 9.0, 7.0, 6.0, 8.0, 9.0, 6.0, 6.0, 1.0, 1.0, 2.0, 8.0, 7.0, 5.0]
global b_x = 5
global d_y = [4.0, 5.0, 4.0, 10.0, 4.0, 2.0, 5.0, 7.0, 1.0, 6.0, 6.0, 10.0, 7.0, 1.0, 4.0, 9.0, 5.0, 9.0, 1.0, 10.0, 5.0, 5.0, 6.0, 6.0, 3.0, 1.0, 10.0, 2.0, 5.0, 5.0, 10.0, 9.0, 3.0, 8.0, 1.0, 7.0, 4.0, 8.0, 6.0, 10.0, 8.0, 9.0, 8.0, 2.0, 6.0, 9.0, 8.0, 9.0, 9.0, 8.0, 3.0, 9.0, 8.0, 3.0, 2.0, 7.0, 5.0, 5.0, 9.0, 1.0, 3.0, 6.0, 10.0, 4.0, 9.0, 8.0, 8.0, 8.0, 7.0, 7.0, 6.0, 7.0, 9.0, 5.0, 3.0, 1.0, 7.0, 9.0, 7.0, 3.0, 5.0, 4.0, 6.0, 4.0, 7.0, 7.0, 10.0, 3.0, 1.0, 8.0, 8.0, 1.0, 9.0, 3.0, 2.0, 8.0, 3.0, 7.0, 10.0, 5.0, 9.0, 10.0, 4.0, 8.0, 1.0, 4.0, 9.0, 7.0, 5.0, 9.0, 8.0, 9.0, 4.0, 4.0, 2.0, 6.0, 2.0, 5.0, 6.0, 1.0, 7.0]
global b_y = 10
global p = [0.769, 0.974, 0.688, 0.221, 0.846, 0.516, 0.68, 0.846, 0.269, 0.466, 0.437, 0.737, 0.031, 0.392, 0.728, 0.677, 0.24, 0.537, 0.513, 0.86, 0.863, 0.832, 0.775, 0.841, 0.71, 0.777, 0.778, 0.456, 0.423, 0.94, 0.065, 0.131, 0.235, 0.329, 0.856, 0.737, 0.952, 0.329, 0.984, 0.078, 0.021, 0.494, 0.554, 0.725, 0.187, 0.102, 0.897, 0.213, 0.695, 0.154, 0.827, 0.308, 0.438, 0.866, 0.875, 0.556, 0.954, 0.755, 0.298, 0.261, 0.807, 0.895, 0.787, 0.681, 0.725, 0.076, 0.871, 0.668, 0.07, 0.939, 0.659, 0.94, 0.801, 0.381, 0.182, 0.729, 0.94, 0.599, 0.533, 0.233, 0.081, 0.96, 0.726, 0.54, 0.203, 0.364, 0.423, 0.042, 0.758, 0.635, 0.369, 0.928, 0.4, 0.94, 0.762, 0.325, 0.844, 0.14, 0.087, 0.649, 0.118, 0.071, 0.815, 0.238, 0.417, 0.187, 0.895, 0.211, 0.826, 0.723, 0.421, 0.238, 0.365, 0.666, 0.803, 0.029, 0.611, 0.5, 0.472, 0.085, 0.651]
global q = [0.876, 0.985, 0.757, 0.476, 0.919, 0.623, 0.834, 0.865, 0.349, 0.716, 0.946, 0.98, 0.84, 0.493, 0.733, 0.852, 0.491, 0.789, 0.842, 0.922, 0.891, 0.973, 0.968, 0.88, 0.913, 0.803, 0.843, 0.726, 0.49, 0.973, 0.475, 0.284, 0.63, 0.695, 0.956, 0.954, 0.987, 0.969, 0.991, 0.099, 0.467, 0.654, 0.998, 0.93, 0.841, 0.698, 0.913, 0.316, 0.94, 0.918, 0.961, 0.539, 0.648, 0.913, 0.914, 0.807, 0.978, 0.924, 0.601, 0.82, 0.991, 0.957, 0.816, 0.799, 0.955, 0.936, 0.92, 0.856, 0.681, 0.961, 0.991, 0.973, 0.818, 0.448, 0.808, 0.888, 0.979, 0.892, 0.828, 0.25, 0.499, 0.964, 0.982, 0.712, 0.537, 0.564, 0.798, 0.481, 0.966, 0.716, 0.615, 0.937, 0.706, 0.993, 0.846, 0.583, 0.973, 0.218, 0.997, 0.747, 0.887, 0.912, 0.989, 0.572, 0.52, 0.609, 0.979, 0.769, 0.855, 0.883, 0.82, 0.969, 0.548, 0.967, 0.972, 0.288, 0.666, 0.511, 0.516, 0.217, 0.945]
global origin = 1
global destination = 35