global arcs = [1 14; 1 17; 1 23; 1 36; 2 10; 3 17; 3 19; 3 22; 3 39; 4 12; 4 15; 4 27; 4 29; 4 38; 5 4; 5 8; 5 13; 5 27; 6 8; 6 13; 6 39; 7 19; 7 24; 7 36; 7 40; 8 19; 8 34; 8 35; 9 7; 9 10; 9 15; 9 16; 9 22; 9 25; 10 4; 10 16; 10 34; 11 7; 11 12; 11 13; 11 17; 11 18; 11 23; 11 28; 12 3; 12 13; 12 20; 12 23; 12 25; 12 30; 12 31; 12 36; 13 9; 13 10; 13 16; 13 20; 13 24; 14 11; 14 16; 14 17; 14 30; 14 35; 15 18; 16 9; 16 24; 17 4; 17 13; 18 17; 18 30; 19 3; 19 9; 19 25; 20 13; 20 21; 20 33; 20 37; 21 6; 21 12; 21 23; 21 33; 21 35; 22 10; 22 16; 22 24; 23 13; 23 30; 23 31; 23 36; 24 4; 24 13; 24 20; 24 32; 25 2; 25 6; 25 18; 25 28; 25 31; 26 11; 26 23; 26 31; 27 9; 27 18; 27 31; 28 3; 28 7; 28 12; 28 16; 28 25; 28 36; 29 24; 29 27; 30 15; 31 4; 31 11; 31 21; 31 24; 31 39; 32 9; 32 35; 32 36; 32 37; 33 3; 33 5; 33 23; 33 28; 34 5; 34 11; 34 26; 34 28; 35 6; 35 11; 35 27; 35 32; 35 37; 36 6; 36 10; 36 21; 36 23; 36 28; 36 31; 36 38; 37 12; 37 13; 37 28; 38 25; 38 26; 38 29; 39 30; 39 40]
global d_x = [6.0, 2.0, 1.0, 1.0, 4.0, 5.0, 2.0, 3.0, 5.0, 3.0, 8.0, 6.0, 8.0, 10.0, 10.0, 8.0, 10.0, 9.0, 5.0, 9.0, 2.0, 2.0, 9.0, 9.0, 4.0, 6.0, 9.0, 5.0, 6.0, 10.0, 3.0, 4.0, 3.0, 2.0, 10.0, 10.0, 4.0, 5.0, 4.0, 2.0, 6.0, 6.0, 9.0, 7.0, 2.0, 5.0, 2.0, 1.0, 1.0, 10.0, 10.0, 3.0, 9.0, 4.0, 5.0, 9.0, 9.0, 10.0, 1.0, 1.0, 6.0, 7.0, 9.0, 9.0, 1.0, 6.0, 1.0, 3.0, 1.0, 6.0, 1.0, 9.0, 9.0, 8.0, 2.0, 1.0, 10.0, 3.0, 6.0, 1.0, 10.0, 6.0, 8.0, 2.0, 6.0, 8.0, 2.0, 9.0, 6.0, 1.0, 6.0, 2.0, 9.0, 2.0, 1.0, 10.0, 5.0, 2.0, 7.0, 2.0, 1.0, 2.0, 6.0, 2.0, 10.0, 3.0, 4.0, 5.0, 7.0, 5.0, 7.0, 8.0, 8.0, 8.0, 1.0, 5.0, 6.0, 1.0, 10.0, 9.0, 3.0, 6.0, 9.0, 3.0, 3.0, 3.0, 8.0, 10.0, 4.0, 9.0, 10.0, 8.0, 8.0, 9.0, 4.0, 2.0, 2.0, 10.0, 1.0, 9.0, 2.0, 9.0, 2.0, 8.0, 3.0, 9.0, 6.0, 9.0, 4.0]
global b_x = 5
global d_y = [3.0, 4.0, 2.0, 6.0, 8.0, 9.0, 3.0, 10.0, 5.0, 9.0, 4.0, 5.0, 8.0, 5.0, 9.0, 8.0, 2.0, 7.0, 3.0, 10.0, 6.0, 8.0, 7.0, 3.0, 2.0, 1.0, 2.0, 8.0, 1.0, 4.0, 5.0, 9.0, 9.0, 8.0, 8.0, 8.0, 1.0, 10.0, 1.0, 9.0, 9.0, 4.0, 5.0, 3.0, 8.0, 9.0, 6.0, 10.0, 8.0, 10.0, 4.0, 10.0, 7.0, 5.0, 5.0, 2.0, 2.0, 5.0, 9.0, 9.0, 5.0, 3.0, 6.0, 10.0, 4.0, 4.0, 3.0, 6.0, 2.0, 7.0, 5.0, 9.0, 3.0, 10.0, 8.0, 7.0, 10.0, 1.0, 4.0, 3.0, 10.0, 2.0, 7.0, 10.0, 7.0, 4.0, 10.0, 9.0, 4.0, 3.0, 8.0, 4.0, 2.0, 1.0, 1.0, 6.0, 10.0, 4.0, 2.0, 8.0, 8.0, 7.0, 9.0, 10.0, 2.0, 3.0, 4.0, 9.0, 5.0, 9.0, 1.0, 3.0, 7.0, 2.0, 8.0, 5.0, 3.0, 1.0, 5.0, 8.0, 3.0, 10.0, 3.0, 6.0, 10.0, 9.0, 2.0, 6.0, 1.0, 9.0, 3.0, 3.0, 5.0, 5.0, 3.0, 6.0, 2.0, 1.0, 4.0, 5.0, 3.0, 2.0, 7.0, 3.0, 1.0, 9.0, 8.0, 5.0, 2.0]
global b_y = 10
global p = [0.043, 0.461, 0.312, 0.362, 0.63, 0.97, 0.356, 0.474, 0.877, 0.028, 0.062, 0.413, 0.446, 0.509, 0.906, 0.805, 0.975, 0.607, 0.649, 0.069, 0.945, 0.56, 0.335, 0.329, 0.529, 0.839, 0.438, 0.999, 0.527, 0.626, 0.782, 0.837, 0.613, 0.425, 0.772, 0.419, 0.766, 0.915, 0.128, 0.83, 0.008, 0.649, 0.403, 0.897, 0.975, 0.776, 0.155, 0.093, 0.762, 0.897, 0.917, 0.951, 0.596, 0.883, 0.323, 0.725, 0.133, 0.504, 0.91, 0.532, 0.394, 0.579, 0.471, 0.012, 0.072, 0.124, 0.754, 0.612, 0.487, 0.373, 0.291, 0.566, 0.836, 0.531, 0.276, 0.849, 0.017, 0.942, 0.416, 0.984, 0.482, 0.953, 0.218, 0.644, 0.521, 0.568, 0.484, 0.279, 0.009, 0.886, 0.444, 0.003, 0.997, 0.199, 0.337, 0.919, 0.265, 0.953, 0.27, 0.06, 0.645, 0.226, 0.198, 0.011, 0.869, 0.65, 0.495, 0.29, 0.194, 0.873, 0.779, 0.026, 0.893, 0.031, 0.885, 0.55, 0.222, 0.033, 0.399, 0.464, 0.598, 0.482, 0.493, 0.529, 0.365, 0.904, 0.819, 0.294, 0.22, 0.523, 0.711, 0.963, 0.078, 0.156, 0.365, 0.305, 0.551, 0.555, 0.959, 0.584, 0.836, 0.352, 0.058, 0.313, 0.727, 0.383, 0.117, 0.614, 0.683]
global q = [0.34, 0.924, 0.62, 0.864, 0.998, 0.988, 0.415, 0.745, 0.905, 0.363, 0.247, 0.434, 0.677, 0.958, 0.953, 0.951, 0.984, 0.679, 0.983, 0.403, 0.973, 0.754, 0.72, 0.437, 0.937, 0.863, 0.615, 0.999, 0.878, 0.64, 0.909, 0.99, 0.919, 0.52, 0.794, 0.425, 0.79, 0.995, 0.692, 0.915, 0.898, 0.876, 0.775, 0.998, 0.977, 0.852, 0.501, 0.616, 0.851, 0.928, 0.925, 0.954, 0.762, 0.949, 0.611, 0.743, 0.17, 0.662, 0.985, 0.915, 0.864, 0.95, 0.563, 0.342, 0.704, 0.404, 0.902, 0.821, 0.839, 0.977, 0.813, 0.689, 0.996, 0.625, 0.929, 0.876, 0.086, 0.991, 0.742, 0.991, 0.934, 0.977, 0.7, 0.89, 0.999, 0.657, 0.487, 0.462, 0.964, 0.887, 0.956, 0.963, 0.997, 0.615, 0.798, 0.993, 0.974, 0.995, 0.491, 0.101, 0.725, 0.293, 0.847, 0.896, 0.916, 0.895, 0.898, 0.618, 0.902, 0.911, 0.861, 0.924, 0.901, 0.416, 0.994, 0.854, 0.335, 0.66, 0.991, 0.726, 0.886, 0.817, 0.543, 0.706, 0.502, 0.95, 0.974, 0.867, 0.886, 0.807, 0.785, 0.966, 0.09, 0.479, 0.851, 0.473, 0.61, 0.613, 0.976, 0.738, 0.927, 0.547, 0.489, 0.876, 0.872, 0.921, 0.136, 0.769, 0.92]
global origin = 1
global destination = 40