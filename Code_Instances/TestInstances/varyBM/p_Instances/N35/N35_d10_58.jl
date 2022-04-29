global arcs = [1 3; 2 10; 2 18; 2 26; 2 29; 2 30; 3 15; 3 16; 3 27; 4 13; 5 10; 5 27; 5 33; 6 12; 6 19; 7 3; 7 4; 7 5; 7 14; 7 33; 8 2; 8 6; 9 2; 9 3; 9 5; 9 13; 9 21; 9 28; 10 6; 10 28; 10 29; 11 6; 11 29; 11 30; 11 32; 12 7; 12 31; 13 7; 13 15; 13 20; 14 4; 14 9; 14 10; 15 18; 16 3; 16 19; 16 35; 17 4; 17 18; 17 24; 17 25; 17 31; 18 2; 18 13; 19 4; 19 21; 19 23; 19 27; 20 4; 20 12; 20 15; 20 18; 20 24; 20 33; 21 11; 21 17; 21 19; 21 23; 21 31; 22 7; 22 14; 22 19; 22 29; 23 4; 23 7; 23 13; 23 20; 23 25; 23 26; 23 29; 23 34; 24 5; 24 8; 24 23; 24 33; 25 4; 25 6; 25 8; 25 22; 25 23; 25 26; 26 8; 26 13; 26 24; 26 34; 27 9; 27 10; 27 16; 27 30; 28 7; 28 8; 28 10; 28 17; 28 22; 29 24; 29 27; 29 30; 30 7; 30 12; 30 13; 30 32; 30 35; 31 9; 31 13; 31 22; 31 25; 31 28; 31 32; 31 34; 32 2; 32 24; 32 35; 33 11; 33 35; 34 12; 34 13; 34 19]
global d_x = [6.0, 4.0, 8.0, 2.0, 8.0, 8.0, 9.0, 1.0, 8.0, 6.0, 1.0, 8.0, 3.0, 1.0, 5.0, 1.0, 2.0, 3.0, 1.0, 6.0, 10.0, 1.0, 1.0, 4.0, 1.0, 5.0, 7.0, 7.0, 9.0, 3.0, 3.0, 5.0, 10.0, 6.0, 7.0, 10.0, 9.0, 8.0, 6.0, 1.0, 1.0, 8.0, 5.0, 7.0, 5.0, 6.0, 1.0, 9.0, 2.0, 9.0, 9.0, 3.0, 1.0, 10.0, 7.0, 9.0, 10.0, 7.0, 6.0, 5.0, 10.0, 7.0, 4.0, 3.0, 1.0, 5.0, 6.0, 3.0, 5.0, 3.0, 9.0, 9.0, 4.0, 7.0, 4.0, 9.0, 9.0, 7.0, 4.0, 8.0, 9.0, 5.0, 3.0, 2.0, 7.0, 5.0, 6.0, 9.0, 1.0, 6.0, 4.0, 9.0, 6.0, 6.0, 10.0, 2.0, 5.0, 6.0, 1.0, 4.0, 10.0, 8.0, 7.0, 7.0, 1.0, 4.0, 3.0, 5.0, 9.0, 5.0, 7.0, 10.0, 1.0, 8.0, 8.0, 9.0, 2.0, 6.0, 4.0, 6.0, 5.0, 2.0, 5.0, 6.0, 6.0, 8.0, 8.0]
global b_x = 5
global d_y = [1.0, 8.0, 9.0, 9.0, 6.0, 9.0, 5.0, 10.0, 6.0, 9.0, 10.0, 5.0, 4.0, 1.0, 9.0, 3.0, 6.0, 10.0, 8.0, 5.0, 7.0, 10.0, 2.0, 7.0, 1.0, 2.0, 4.0, 3.0, 6.0, 4.0, 6.0, 8.0, 8.0, 10.0, 3.0, 10.0, 4.0, 7.0, 5.0, 6.0, 4.0, 9.0, 4.0, 8.0, 9.0, 3.0, 1.0, 2.0, 5.0, 3.0, 6.0, 5.0, 9.0, 2.0, 8.0, 6.0, 3.0, 7.0, 7.0, 3.0, 5.0, 3.0, 6.0, 3.0, 1.0, 7.0, 9.0, 3.0, 1.0, 1.0, 4.0, 1.0, 8.0, 4.0, 9.0, 9.0, 6.0, 9.0, 8.0, 7.0, 5.0, 9.0, 4.0, 8.0, 10.0, 2.0, 10.0, 6.0, 1.0, 4.0, 10.0, 9.0, 2.0, 8.0, 6.0, 4.0, 4.0, 6.0, 10.0, 3.0, 7.0, 5.0, 10.0, 10.0, 7.0, 9.0, 2.0, 10.0, 8.0, 4.0, 5.0, 2.0, 5.0, 3.0, 9.0, 4.0, 3.0, 2.0, 7.0, 8.0, 6.0, 4.0, 4.0, 10.0, 7.0, 5.0, 4.0]
global b_y = 10
global p = [0.791, 0.884, 0.89, 0.534, 0.459, 0.646, 0.656, 0.17, 0.02, 0.765, 0.299, 0.439, 0.224, 0.148, 0.056, 0.34, 0.504, 0.054, 0.544, 0.816, 0.087, 0.951, 0.606, 0.831, 0.948, 0.152, 0.484, 0.694, 0.678, 0.899, 0.213, 0.435, 0.733, 0.468, 0.889, 0.171, 0.67, 0.66, 0.763, 0.587, 0.28, 0.562, 0.46, 0.694, 0.294, 0.316, 0.385, 0.505, 0.895, 0.758, 0.496, 0.142, 0.313, 0.506, 0.805, 0.312, 0.06, 0.195, 0.53, 0.97, 0.926, 0.213, 0.022, 0.103, 0.859, 0.647, 0.921, 0.689, 0.314, 0.936, 0.007, 0.39, 0.761, 0.566, 0.352, 0.763, 0.543, 0.516, 0.35, 0.162, 0.12, 0.649, 0.894, 0.297, 0.139, 0.78, 0.876, 0.317, 0.375, 0.276, 0.917, 0.94, 0.84, 0.181, 0.04, 0.482, 0.52, 0.263, 0.251, 0.594, 0.113, 0.883, 0.871, 0.522, 0.671, 0.644, 0.543, 0.346, 0.305, 0.537, 0.175, 0.52, 0.612, 0.619, 0.53, 0.059, 0.263, 0.286, 0.635, 0.005, 0.602, 0.311, 0.112, 0.22, 0.053, 0.467, 0.282]
global q = [0.821, 0.915, 0.924, 0.673, 0.602, 0.893, 0.926, 0.474, 0.698, 0.774, 0.324, 0.595, 0.685, 0.873, 0.813, 0.707, 0.811, 0.468, 0.789, 0.902, 0.997, 0.976, 0.628, 0.867, 0.948, 0.404, 0.51, 0.811, 0.8, 0.969, 0.885, 0.561, 0.961, 0.988, 0.98, 0.969, 0.989, 0.742, 0.803, 0.822, 0.716, 0.803, 0.886, 0.925, 0.484, 0.472, 0.62, 0.568, 0.943, 0.863, 0.936, 0.296, 0.811, 0.542, 0.879, 0.448, 0.769, 0.33, 0.962, 0.98, 0.973, 0.509, 0.854, 0.944, 0.869, 0.916, 0.954, 0.891, 0.718, 0.985, 0.688, 0.773, 0.769, 0.939, 0.4, 0.954, 0.758, 0.917, 0.429, 0.612, 0.136, 0.874, 0.961, 0.849, 0.428, 0.795, 0.921, 0.346, 0.752, 0.461, 0.998, 0.951, 0.86, 0.363, 0.989, 0.671, 0.766, 0.562, 0.571, 0.905, 0.873, 0.885, 0.923, 0.588, 0.689, 0.756, 0.782, 0.358, 0.707, 0.617, 0.246, 0.739, 0.805, 0.789, 0.627, 0.226, 0.697, 0.962, 0.766, 0.092, 0.952, 0.738, 0.591, 0.383, 0.719, 0.938, 0.36]
global origin = 1
global destination = 35
