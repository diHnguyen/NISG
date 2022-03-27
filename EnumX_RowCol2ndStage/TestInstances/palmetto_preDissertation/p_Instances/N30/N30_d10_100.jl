global arcs = [1 10; 1 26; 2 23; 2 28; 3 9; 3 13; 3 20; 4 10; 4 19; 4 23; 4 27; 5 23; 5 24; 5 28; 6 2; 6 9; 6 14; 6 15; 6 26; 6 30; 7 2; 7 3; 7 4; 7 12; 8 3; 8 22; 8 29; 9 12; 9 15; 10 24; 10 25; 11 3; 11 6; 11 15; 11 24; 12 9; 12 17; 12 28; 13 27; 14 8; 14 24; 15 4; 15 17; 15 18; 15 24; 15 30; 16 2; 16 4; 16 6; 16 8; 16 17; 16 18; 16 22; 17 3; 17 10; 17 13; 17 15; 17 20; 18 4; 18 5; 18 17; 18 19; 18 20; 18 30; 19 13; 19 14; 19 24; 20 3; 20 21; 21 2; 21 20; 22 8; 22 10; 22 12; 22 28; 23 7; 23 11; 23 22; 24 16; 25 4; 25 5; 25 9; 26 4; 26 19; 26 20; 26 27; 26 30; 27 6; 27 16; 28 3; 28 27; 29 8; 29 11; 29 13; 29 16]
global d_x = [9.0, 10.0, 6.0, 3.0, 4.0, 9.0, 1.0, 5.0, 2.0, 7.0, 6.0, 9.0, 2.0, 6.0, 8.0, 7.0, 4.0, 1.0, 1.0, 7.0, 5.0, 2.0, 2.0, 2.0, 5.0, 1.0, 5.0, 9.0, 6.0, 2.0, 6.0, 9.0, 2.0, 6.0, 10.0, 4.0, 9.0, 7.0, 3.0, 9.0, 4.0, 3.0, 4.0, 7.0, 8.0, 2.0, 1.0, 10.0, 9.0, 2.0, 10.0, 1.0, 8.0, 9.0, 8.0, 4.0, 9.0, 9.0, 2.0, 3.0, 9.0, 7.0, 7.0, 3.0, 5.0, 5.0, 7.0, 2.0, 3.0, 7.0, 8.0, 5.0, 2.0, 5.0, 5.0, 2.0, 6.0, 4.0, 9.0, 4.0, 10.0, 10.0, 5.0, 7.0, 5.0, 10.0, 9.0, 2.0, 9.0, 10.0, 9.0, 8.0, 8.0, 9.0, 4.0]
global b_x = 5
global d_y = [1.0, 4.0, 1.0, 1.0, 9.0, 3.0, 9.0, 10.0, 8.0, 7.0, 1.0, 3.0, 6.0, 6.0, 6.0, 10.0, 9.0, 2.0, 10.0, 6.0, 2.0, 5.0, 2.0, 6.0, 7.0, 2.0, 4.0, 2.0, 8.0, 4.0, 1.0, 3.0, 1.0, 3.0, 2.0, 2.0, 1.0, 8.0, 3.0, 7.0, 2.0, 6.0, 4.0, 10.0, 3.0, 3.0, 10.0, 7.0, 2.0, 5.0, 6.0, 10.0, 4.0, 3.0, 2.0, 7.0, 10.0, 3.0, 3.0, 8.0, 6.0, 8.0, 8.0, 3.0, 8.0, 3.0, 9.0, 9.0, 7.0, 10.0, 3.0, 9.0, 2.0, 8.0, 10.0, 4.0, 2.0, 10.0, 1.0, 10.0, 2.0, 4.0, 10.0, 3.0, 9.0, 5.0, 9.0, 1.0, 7.0, 10.0, 10.0, 6.0, 5.0, 6.0, 7.0]
global b_y = 10
global p = [0.212, 0.985, 0.985, 0.272, 0.669, 0.363, 0.024, 0.255, 0.309, 0.817, 0.482, 0.45, 0.334, 0.602, 0.116, 0.44, 0.059, 0.076, 0.671, 0.137, 0.631, 0.533, 0.714, 0.51, 0.071, 0.772, 0.603, 0.792, 0.857, 0.916, 0.563, 0.213, 0.525, 0.393, 0.559, 0.309, 0.21, 0.833, 0.267, 0.279, 0.572, 0.568, 0.582, 0.037, 0.977, 0.085, 0.677, 0.946, 0.04, 0.558, 0.581, 0.12, 0.272, 0.674, 0.558, 0.347, 0.83, 0.134, 0.042, 0.803, 0.45, 0.567, 0.172, 0.232, 0.18, 0.777, 0.261, 0.957, 0.342, 0.941, 0.509, 0.267, 0.477, 0.626, 0.472, 0.915, 0.567, 0.539, 0.242, 0.608, 0.353, 0.071, 0.771, 0.256, 0.039, 0.527, 0.943, 0.66, 0.162, 0.171, 0.866, 0.343, 0.139, 0.08, 0.592]
global q = [0.362, 0.991, 0.988, 0.315, 0.949, 0.91, 0.439, 0.631, 0.419, 0.968, 0.868, 0.99, 0.779, 0.752, 0.701, 0.57, 0.936, 0.976, 0.763, 0.196, 0.8, 0.692, 0.83, 0.512, 0.814, 0.87, 0.968, 0.873, 0.946, 0.931, 0.919, 0.606, 0.904, 0.954, 0.671, 0.577, 0.523, 0.858, 0.687, 0.318, 0.883, 0.571, 0.852, 0.914, 0.982, 0.392, 0.944, 0.967, 0.781, 0.764, 0.659, 0.581, 0.945, 0.771, 0.977, 0.999, 0.909, 0.15, 0.41, 0.883, 0.645, 0.925, 0.874, 0.384, 0.674, 0.973, 0.326, 0.957, 0.353, 0.951, 0.824, 0.594, 0.962, 0.672, 0.828, 0.938, 0.736, 0.916, 0.351, 0.866, 0.72, 0.971, 0.969, 0.95, 0.473, 0.57, 0.975, 0.707, 0.402, 0.914, 0.983, 0.609, 0.594, 0.961, 0.714]
global origin = 1
global destination = 30