global arcs = [1 4; 1 19; 1 23; 1 31; 1 40; 2 4; 2 5; 2 7; 2 16; 2 34; 2 39; 3 4; 3 14; 3 17; 3 18; 3 32; 3 34; 4 10; 4 38; 5 6; 5 10; 5 12; 5 24; 5 29; 5 30; 5 38; 6 5; 6 8; 6 24; 6 26; 6 36; 7 13; 7 15; 7 16; 7 20; 7 23; 7 28; 7 29; 8 13; 8 16; 8 22; 8 34; 8 39; 9 2; 9 12; 9 17; 9 18; 9 22; 10 9; 10 24; 11 4; 11 6; 11 9; 11 16; 11 17; 11 26; 11 33; 11 38; 12 16; 12 21; 12 24; 12 38; 13 4; 13 5; 13 8; 13 18; 13 25; 13 32; 13 35; 13 40; 14 5; 14 8; 14 27; 14 36; 15 9; 15 16; 15 26; 15 36; 16 7; 16 21; 16 22; 16 24; 16 28; 17 16; 17 24; 17 30; 17 33; 18 8; 18 11; 18 24; 18 38; 18 40; 19 23; 19 36; 20 4; 20 8; 20 23; 20 29; 20 39; 21 29; 21 32; 21 33; 21 38; 22 7; 22 9; 22 15; 22 27; 22 33; 23 2; 23 3; 23 31; 24 7; 24 11; 24 20; 24 38; 25 4; 25 7; 25 10; 25 21; 25 24; 25 39; 26 5; 26 18; 26 25; 26 38; 27 4; 27 6; 27 7; 27 10; 27 11; 27 37; 28 9; 28 12; 28 15; 29 16; 29 25; 29 28; 29 32; 29 33; 29 35; 30 4; 30 23; 30 26; 30 35; 30 36; 31 14; 31 15; 31 24; 32 6; 32 23; 32 29; 32 33; 32 35; 32 38; 33 34; 34 21; 34 27; 35 10; 36 13; 37 27; 37 30; 38 5; 38 6; 38 13; 38 18; 38 25; 38 27; 38 34; 39 30]
global d_x = [4.0, 5.0, 7.0, 4.0, 6.0, 9.0, 2.0, 8.0, 8.0, 2.0, 5.0, 4.0, 8.0, 2.0, 2.0, 1.0, 7.0, 6.0, 2.0, 5.0, 3.0, 2.0, 5.0, 3.0, 8.0, 10.0, 8.0, 2.0, 8.0, 5.0, 8.0, 4.0, 4.0, 5.0, 1.0, 2.0, 6.0, 8.0, 8.0, 9.0, 9.0, 3.0, 7.0, 4.0, 3.0, 8.0, 9.0, 10.0, 10.0, 2.0, 8.0, 8.0, 10.0, 9.0, 7.0, 9.0, 2.0, 5.0, 10.0, 10.0, 1.0, 4.0, 8.0, 9.0, 10.0, 8.0, 10.0, 7.0, 7.0, 5.0, 7.0, 9.0, 3.0, 4.0, 4.0, 5.0, 6.0, 9.0, 2.0, 6.0, 4.0, 2.0, 6.0, 4.0, 10.0, 4.0, 2.0, 4.0, 7.0, 5.0, 1.0, 8.0, 4.0, 1.0, 9.0, 8.0, 9.0, 8.0, 1.0, 5.0, 9.0, 6.0, 2.0, 7.0, 5.0, 1.0, 8.0, 3.0, 4.0, 7.0, 6.0, 6.0, 10.0, 2.0, 9.0, 1.0, 6.0, 3.0, 8.0, 4.0, 4.0, 5.0, 3.0, 7.0, 9.0, 9.0, 7.0, 4.0, 7.0, 2.0, 9.0, 4.0, 4.0, 7.0, 8.0, 9.0, 1.0, 3.0, 1.0, 10.0, 1.0, 9.0, 1.0, 2.0, 9.0, 7.0, 3.0, 2.0, 4.0, 9.0, 4.0, 2.0, 1.0, 10.0, 8.0, 5.0, 4.0, 5.0, 1.0, 3.0, 6.0, 4.0, 8.0, 7.0, 4.0, 10.0, 6.0, 8.0, 9.0]
global b_x = 5
global d_y = [3.0, 3.0, 9.0, 3.0, 1.0, 4.0, 2.0, 6.0, 2.0, 4.0, 2.0, 1.0, 7.0, 6.0, 1.0, 5.0, 5.0, 4.0, 5.0, 3.0, 8.0, 7.0, 7.0, 8.0, 8.0, 3.0, 7.0, 6.0, 9.0, 6.0, 10.0, 5.0, 4.0, 7.0, 10.0, 1.0, 3.0, 7.0, 8.0, 6.0, 6.0, 10.0, 6.0, 9.0, 10.0, 8.0, 6.0, 2.0, 4.0, 7.0, 1.0, 3.0, 4.0, 7.0, 10.0, 1.0, 6.0, 5.0, 9.0, 5.0, 8.0, 8.0, 4.0, 2.0, 5.0, 8.0, 3.0, 2.0, 10.0, 3.0, 7.0, 8.0, 2.0, 3.0, 9.0, 7.0, 10.0, 4.0, 10.0, 4.0, 3.0, 6.0, 7.0, 1.0, 6.0, 4.0, 5.0, 3.0, 2.0, 7.0, 6.0, 6.0, 7.0, 4.0, 6.0, 10.0, 2.0, 1.0, 5.0, 3.0, 7.0, 6.0, 7.0, 8.0, 6.0, 4.0, 2.0, 8.0, 7.0, 2.0, 9.0, 7.0, 1.0, 10.0, 5.0, 5.0, 7.0, 10.0, 10.0, 1.0, 9.0, 2.0, 1.0, 1.0, 5.0, 3.0, 3.0, 1.0, 7.0, 3.0, 9.0, 4.0, 3.0, 5.0, 2.0, 9.0, 9.0, 2.0, 1.0, 1.0, 2.0, 9.0, 3.0, 4.0, 2.0, 6.0, 5.0, 10.0, 7.0, 10.0, 5.0, 8.0, 6.0, 4.0, 1.0, 8.0, 3.0, 9.0, 3.0, 1.0, 2.0, 8.0, 6.0, 9.0, 10.0, 9.0, 4.0, 7.0, 2.0]
global b_y = 10
global p = [0.419, 0.239, 0.731, 0.119, 0.884, 0.767, 0.797, 0.125, 0.734, 0.556, 0.092, 0.776, 0.355, 0.705, 0.638, 0.589, 0.922, 0.62, 0.144, 0.453, 0.988, 0.157, 0.838, 0.934, 0.936, 0.945, 0.684, 0.844, 0.496, 0.956, 0.077, 0.71, 0.959, 0.992, 0.601, 0.583, 0.256, 0.966, 0.352, 0.747, 0.42, 0.544, 0.296, 0.501, 0.534, 0.587, 0.943, 0.946, 0.288, 0.331, 0.123, 0.342, 0.405, 0.915, 0.895, 0.754, 0.455, 0.55, 0.198, 0.204, 0.382, 0.896, 0.305, 0.986, 0.602, 0.133, 0.646, 0.517, 0.694, 0.406, 0.082, 0.279, 0.945, 0.996, 0.398, 0.901, 0.073, 0.394, 0.086, 0.944, 0.248, 0.773, 0.201, 0.806, 0.489, 0.757, 0.298, 0.141, 0.121, 0.543, 0.736, 0.007, 0.987, 0.624, 0.893, 0.904, 0.57, 0.303, 0.507, 0.661, 0.289, 0.966, 0.587, 0.692, 0.455, 0.234, 0.35, 0.943, 0.456, 0.227, 0.963, 0.054, 0.161, 0.369, 0.759, 0.119, 0.496, 0.364, 0.08, 0.543, 0.539, 0.197, 0.137, 0.56, 0.964, 0.22, 0.336, 0.806, 0.967, 0.903, 0.631, 0.355, 0.081, 0.348, 0.639, 0.917, 0.653, 0.236, 0.596, 0.47, 0.667, 0.779, 0.766, 0.472, 0.847, 0.913, 0.404, 0.094, 0.353, 0.776, 0.753, 0.191, 0.842, 0.248, 0.73, 0.357, 0.401, 0.464, 0.977, 0.429, 0.328, 0.905, 0.893, 0.822, 0.738, 0.414, 0.828, 0.895, 0.454]
global q = [0.556, 0.309, 0.792, 0.527, 0.962, 0.926, 0.967, 0.507, 0.914, 0.669, 0.781, 0.832, 0.942, 0.997, 0.782, 0.997, 0.959, 0.892, 0.574, 0.692, 0.999, 0.886, 0.846, 0.979, 0.963, 0.961, 0.918, 0.98, 0.988, 0.964, 0.127, 0.789, 0.968, 0.999, 0.708, 0.603, 0.761, 0.967, 0.695, 0.755, 0.998, 0.901, 0.95, 0.723, 0.805, 0.781, 0.959, 0.994, 0.467, 0.405, 0.823, 0.928, 0.501, 0.945, 0.934, 0.921, 0.641, 0.995, 0.602, 0.479, 0.415, 0.897, 0.883, 0.997, 0.866, 0.393, 0.809, 0.916, 0.979, 0.818, 0.216, 0.914, 0.986, 0.998, 0.699, 0.931, 0.599, 0.883, 0.618, 0.996, 0.75, 0.866, 0.935, 0.943, 0.734, 0.927, 0.602, 0.801, 0.368, 0.73, 0.921, 0.469, 0.993, 0.992, 0.97, 0.971, 0.657, 0.644, 0.726, 0.662, 0.457, 0.979, 0.619, 0.962, 0.642, 0.769, 0.61, 0.984, 0.647, 0.443, 0.985, 0.241, 0.448, 0.98, 0.776, 0.164, 0.955, 0.704, 0.872, 0.702, 0.566, 0.756, 0.59, 0.815, 0.971, 0.293, 0.928, 0.88, 0.992, 0.918, 0.917, 0.51, 0.509, 0.643, 0.886, 0.985, 0.661, 0.924, 0.683, 0.779, 0.756, 0.964, 0.775, 0.87, 0.891, 0.972, 0.597, 0.43, 0.786, 0.897, 0.946, 0.784, 0.93, 0.933, 0.95, 0.419, 0.755, 0.959, 0.997, 0.741, 0.347, 0.937, 0.897, 0.885, 0.852, 0.57, 0.938, 0.996, 0.541]
global origin = 1
global destination = 40