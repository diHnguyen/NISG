global arcs = [1 11; 1 25; 1 28; 2 34; 3 10; 3 17; 3 35; 4 10; 4 13; 5 11; 5 12; 5 20; 5 24; 5 34; 5 35; 6 4; 6 8; 6 25; 6 33; 6 35; 7 3; 7 33; 7 37; 8 2; 8 11; 8 12; 8 18; 8 32; 9 8; 9 10; 9 14; 9 17; 9 20; 10 27; 11 4; 11 15; 11 25; 12 13; 12 33; 13 2; 13 5; 13 20; 13 38; 14 15; 14 18; 14 24; 14 29; 14 31; 14 40; 15 27; 15 32; 16 4; 16 18; 16 34; 17 18; 17 19; 17 27; 17 32; 17 34; 18 13; 18 19; 18 28; 18 31; 19 4; 19 12; 19 24; 20 8; 20 10; 20 39; 21 22; 21 24; 21 30; 22 14; 22 19; 22 33; 22 34; 23 7; 23 31; 23 32; 23 40; 24 14; 24 17; 24 29; 25 5; 25 17; 25 33; 25 35; 26 6; 26 21; 26 24; 26 25; 26 29; 26 35; 27 12; 27 36; 28 17; 28 19; 28 25; 28 35; 28 36; 28 37; 29 9; 29 22; 29 40; 30 12; 30 15; 30 35; 30 39; 31 4; 31 13; 31 17; 32 7; 32 10; 32 14; 32 20; 32 27; 32 28; 32 34; 32 35; 32 38; 33 9; 33 16; 33 36; 34 4; 34 20; 34 21; 34 33; 34 40; 35 6; 35 28; 35 40; 36 21; 36 22; 36 37; 36 39; 37 19; 37 26; 38 4; 38 7; 38 14; 38 23; 38 24; 39 3; 39 5; 39 9; 39 22; 39 25; 39 28; 39 31]
global d_x = [8.0, 6.0, 4.0, 4.0, 3.0, 7.0, 4.0, 2.0, 1.0, 1.0, 6.0, 2.0, 6.0, 2.0, 8.0, 10.0, 4.0, 10.0, 10.0, 8.0, 8.0, 2.0, 5.0, 4.0, 10.0, 4.0, 7.0, 6.0, 9.0, 9.0, 7.0, 3.0, 9.0, 7.0, 1.0, 1.0, 3.0, 9.0, 3.0, 3.0, 7.0, 7.0, 6.0, 6.0, 2.0, 4.0, 1.0, 9.0, 9.0, 4.0, 6.0, 7.0, 8.0, 8.0, 10.0, 7.0, 6.0, 2.0, 7.0, 6.0, 7.0, 7.0, 8.0, 6.0, 10.0, 8.0, 6.0, 10.0, 5.0, 2.0, 1.0, 9.0, 7.0, 1.0, 2.0, 9.0, 2.0, 7.0, 3.0, 4.0, 9.0, 5.0, 2.0, 9.0, 4.0, 3.0, 4.0, 2.0, 3.0, 10.0, 6.0, 4.0, 5.0, 3.0, 5.0, 1.0, 10.0, 4.0, 4.0, 1.0, 4.0, 6.0, 4.0, 5.0, 4.0, 5.0, 6.0, 8.0, 4.0, 3.0, 9.0, 7.0, 4.0, 4.0, 10.0, 9.0, 9.0, 3.0, 4.0, 10.0, 3.0, 3.0, 4.0, 4.0, 4.0, 5.0, 1.0, 7.0, 9.0, 3.0, 7.0, 5.0, 3.0, 10.0, 6.0, 2.0, 8.0, 10.0, 3.0, 3.0, 5.0, 4.0, 1.0, 10.0, 6.0, 6.0, 5.0, 7.0, 1.0]
global b_x = 5
global d_y = [5.0, 5.0, 7.0, 3.0, 5.0, 3.0, 1.0, 5.0, 4.0, 10.0, 5.0, 9.0, 3.0, 4.0, 5.0, 10.0, 1.0, 4.0, 3.0, 8.0, 3.0, 3.0, 6.0, 4.0, 8.0, 6.0, 5.0, 8.0, 2.0, 8.0, 8.0, 5.0, 1.0, 10.0, 10.0, 4.0, 10.0, 6.0, 5.0, 6.0, 3.0, 2.0, 2.0, 2.0, 3.0, 4.0, 8.0, 3.0, 6.0, 3.0, 7.0, 7.0, 5.0, 5.0, 4.0, 6.0, 10.0, 5.0, 1.0, 10.0, 7.0, 2.0, 1.0, 7.0, 4.0, 9.0, 6.0, 9.0, 6.0, 10.0, 7.0, 10.0, 5.0, 4.0, 9.0, 3.0, 2.0, 4.0, 2.0, 5.0, 10.0, 6.0, 10.0, 8.0, 1.0, 3.0, 6.0, 3.0, 8.0, 6.0, 4.0, 6.0, 9.0, 3.0, 9.0, 4.0, 6.0, 6.0, 2.0, 7.0, 8.0, 7.0, 10.0, 1.0, 6.0, 1.0, 8.0, 5.0, 4.0, 9.0, 5.0, 5.0, 5.0, 9.0, 10.0, 5.0, 7.0, 2.0, 5.0, 10.0, 5.0, 3.0, 10.0, 10.0, 1.0, 2.0, 1.0, 7.0, 10.0, 7.0, 8.0, 10.0, 1.0, 5.0, 5.0, 6.0, 9.0, 4.0, 9.0, 8.0, 9.0, 4.0, 4.0, 2.0, 9.0, 5.0, 3.0, 6.0, 1.0]
global b_y = 10
global p = [0.451, 0.429, 0.15, 0.536, 0.242, 0.716, 0.468, 0.056, 0.961, 0.213, 0.461, 0.687, 0.299, 0.895, 0.987, 0.961, 0.598, 0.018, 0.683, 0.007, 0.416, 0.756, 0.757, 0.168, 0.981, 0.789, 0.091, 0.477, 0.916, 0.92, 0.531, 0.963, 0.645, 0.799, 0.531, 0.474, 0.369, 0.848, 0.82, 0.829, 0.702, 0.129, 0.488, 0.212, 0.196, 0.766, 0.194, 0.957, 0.367, 0.916, 0.566, 0.016, 0.252, 0.906, 0.652, 0.939, 0.206, 0.756, 0.221, 0.98, 0.564, 0.642, 0.83, 0.762, 0.503, 0.945, 0.096, 0.102, 0.278, 0.578, 0.897, 0.241, 0.892, 0.589, 0.208, 0.91, 0.9, 0.443, 0.43, 0.529, 0.494, 0.317, 0.111, 0.747, 0.23, 0.435, 0.574, 0.947, 0.723, 0.371, 0.16, 0.488, 0.353, 0.461, 0.89, 0.606, 0.479, 0.98, 0.574, 0.342, 0.835, 0.579, 0.623, 0.224, 0.51, 0.805, 0.391, 0.479, 0.532, 0.845, 0.149, 0.664, 0.321, 0.116, 0.536, 0.822, 0.226, 0.762, 0.466, 0.518, 0.836, 0.486, 0.657, 0.64, 0.14, 0.241, 0.256, 0.503, 0.734, 0.647, 0.712, 0.755, 0.733, 0.847, 0.06, 0.684, 0.95, 0.234, 0.683, 0.232, 0.982, 0.703, 0.639, 0.875, 0.28, 0.84, 0.927, 0.39, 0.842]
global q = [0.481, 0.832, 0.667, 0.699, 0.782, 0.756, 0.75, 0.602, 0.962, 0.271, 0.666, 0.721, 0.697, 0.983, 0.995, 0.964, 0.626, 0.432, 0.966, 0.134, 0.828, 0.774, 0.822, 0.376, 0.995, 0.895, 0.455, 0.636, 0.987, 0.932, 0.836, 0.997, 0.672, 0.9, 0.85, 0.632, 0.902, 0.854, 0.84, 0.903, 0.958, 0.2, 0.734, 0.651, 0.649, 0.942, 0.817, 0.967, 0.916, 0.997, 0.833, 0.384, 0.321, 0.942, 0.775, 0.967, 0.507, 0.919, 0.679, 0.983, 0.876, 0.816, 0.892, 0.878, 0.729, 0.982, 0.955, 0.824, 0.836, 0.809, 0.937, 0.702, 0.997, 0.808, 0.501, 0.942, 0.91, 0.994, 0.447, 0.582, 0.598, 0.742, 0.908, 0.954, 0.747, 0.777, 0.663, 0.982, 0.946, 0.769, 0.661, 0.549, 0.573, 0.976, 0.998, 0.885, 0.652, 0.995, 0.67, 0.504, 0.944, 0.706, 0.731, 0.6, 0.852, 0.988, 0.473, 0.637, 0.618, 0.992, 0.193, 0.819, 0.353, 0.3, 0.685, 0.855, 0.579, 0.882, 0.618, 0.796, 0.964, 0.644, 0.924, 0.651, 0.207, 0.551, 0.591, 0.994, 0.843, 0.891, 0.717, 0.975, 0.959, 0.984, 0.858, 0.819, 0.99, 0.267, 0.902, 0.736, 0.989, 0.849, 0.968, 0.996, 0.405, 0.974, 0.952, 0.44, 0.991]
global origin = 1
global destination = 40