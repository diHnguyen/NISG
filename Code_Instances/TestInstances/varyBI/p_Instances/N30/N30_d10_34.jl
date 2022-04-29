global arcs = [1 3; 1 17; 2 26; 2 28; 2 30; 3 13; 3 16; 3 18; 4 5; 4 21; 5 15; 5 26; 6 13; 6 21; 6 27; 7 9; 7 15; 7 27; 7 29; 8 16; 8 24; 8 26; 9 3; 9 10; 9 12; 9 15; 9 18; 10 26; 11 15; 11 30; 12 24; 12 28; 13 2; 13 10; 13 24; 13 28; 14 10; 15 20; 16 15; 16 18; 16 24; 16 26; 17 12; 17 18; 17 19; 18 4; 18 20; 18 22; 19 16; 19 17; 20 8; 20 29; 21 5; 21 14; 21 23; 22 3; 22 11; 22 16; 22 30; 23 3; 23 5; 23 18; 23 25; 24 9; 24 16; 24 17; 24 18; 25 7; 25 26; 26 9; 27 7; 27 15; 28 13; 28 21; 28 22; 28 30; 29 6; 29 12; 29 16; 29 23; 29 25]
global d_x = [5.0, 9.0, 4.0, 2.0, 3.0, 6.0, 6.0, 10.0, 2.0, 6.0, 4.0, 7.0, 1.0, 3.0, 9.0, 6.0, 3.0, 1.0, 6.0, 7.0, 7.0, 7.0, 2.0, 3.0, 5.0, 4.0, 9.0, 4.0, 1.0, 1.0, 5.0, 8.0, 2.0, 8.0, 6.0, 6.0, 10.0, 9.0, 9.0, 8.0, 2.0, 8.0, 1.0, 8.0, 5.0, 2.0, 6.0, 2.0, 1.0, 8.0, 6.0, 6.0, 5.0, 5.0, 1.0, 4.0, 5.0, 10.0, 5.0, 4.0, 8.0, 7.0, 10.0, 5.0, 4.0, 8.0, 9.0, 5.0, 1.0, 1.0, 3.0, 3.0, 4.0, 9.0, 6.0, 9.0, 7.0, 6.0, 10.0, 5.0, 7.0]
global b_x = 5
global d_y = [7.0, 8.0, 1.0, 3.0, 4.0, 1.0, 10.0, 4.0, 4.0, 7.0, 2.0, 7.0, 9.0, 3.0, 10.0, 10.0, 5.0, 1.0, 7.0, 1.0, 4.0, 5.0, 2.0, 2.0, 4.0, 9.0, 1.0, 8.0, 8.0, 3.0, 2.0, 6.0, 2.0, 4.0, 9.0, 3.0, 8.0, 5.0, 9.0, 1.0, 7.0, 8.0, 6.0, 8.0, 10.0, 9.0, 10.0, 7.0, 1.0, 2.0, 1.0, 4.0, 6.0, 10.0, 5.0, 8.0, 10.0, 4.0, 10.0, 4.0, 1.0, 3.0, 3.0, 1.0, 5.0, 7.0, 1.0, 4.0, 6.0, 1.0, 6.0, 7.0, 9.0, 10.0, 2.0, 7.0, 8.0, 2.0, 3.0, 7.0, 1.0]
global b_y = 10
global p = [0.497, 0.232, 0.413, 0.667, 0.979, 0.308, 0.05, 0.12, 0.363, 0.759, 0.421, 0.826, 0.462, 0.403, 0.878, 0.428, 0.145, 0.039, 0.025, 0.069, 0.149, 0.013, 0.235, 0.88, 0.948, 0.073, 0.604, 0.315, 0.905, 0.248, 0.698, 0.327, 0.181, 0.373, 0.705, 0.442, 0.478, 0.906, 0.132, 0.814, 0.042, 0.717, 0.089, 0.648, 0.261, 0.836, 0.322, 0.409, 0.794, 0.791, 0.113, 0.267, 0.859, 0.074, 0.44, 0.071, 0.979, 0.633, 0.536, 0.139, 0.753, 0.618, 0.43, 0.423, 0.782, 0.358, 0.926, 0.773, 0.154, 0.207, 0.186, 0.265, 0.009, 0.22, 0.012, 0.139, 0.711, 0.068, 0.11, 0.188, 0.104]
global q = [0.589, 0.759, 0.798, 0.754, 0.994, 0.795, 0.924, 0.268, 0.497, 0.881, 0.691, 0.871, 0.64, 0.473, 0.888, 0.494, 0.511, 0.92, 0.792, 0.596, 0.679, 0.402, 0.468, 0.906, 0.951, 0.551, 0.717, 0.938, 0.933, 0.776, 0.736, 0.637, 0.715, 0.538, 0.843, 0.733, 0.839, 0.963, 0.499, 0.837, 0.156, 0.951, 0.317, 0.93, 0.928, 0.867, 0.842, 0.585, 0.956, 0.835, 0.325, 0.919, 0.881, 0.248, 0.535, 0.417, 0.983, 0.989, 0.838, 0.443, 0.79, 0.9, 0.617, 0.99, 0.974, 0.726, 0.99, 0.988, 0.586, 0.453, 0.674, 0.413, 0.078, 0.262, 0.642, 0.319, 0.73, 0.239, 0.415, 0.737, 0.678]
global origin = 1
global destination = 30
