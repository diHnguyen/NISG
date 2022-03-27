global arcs = [1 4; 1 17; 1 29; 1 34; 2 16; 2 18; 2 20; 2 35; 2 37; 3 15; 3 21; 3 35; 4 20; 4 37; 5 25; 6 9; 6 14; 6 40; 7 14; 7 16; 7 18; 7 24; 7 40; 8 17; 8 37; 9 2; 9 10; 9 13; 9 21; 9 40; 10 14; 10 17; 10 18; 10 27; 11 24; 11 32; 12 7; 12 17; 13 7; 13 10; 13 11; 13 12; 13 14; 13 36; 14 6; 14 24; 15 23; 15 27; 16 25; 16 26; 17 4; 17 6; 17 32; 18 7; 18 8; 18 10; 18 16; 18 26; 18 34; 19 16; 19 17; 19 30; 20 7; 20 11; 20 32; 21 3; 21 11; 21 18; 21 23; 21 27; 21 32; 21 38; 22 9; 22 10; 22 28; 23 12; 23 33; 23 34; 24 5; 24 11; 24 21; 24 22; 24 26; 24 33; 25 6; 25 26; 25 27; 25 34; 26 2; 26 14; 26 21; 27 11; 27 16; 27 29; 27 31; 27 34; 27 36; 27 37; 28 2; 28 5; 28 27; 29 30; 29 38; 30 3; 30 11; 30 13; 30 17; 30 21; 30 24; 30 35; 30 39; 31 15; 31 37; 32 4; 32 9; 32 10; 32 39; 33 7; 33 14; 34 3; 34 5; 34 20; 34 29; 34 37; 35 19; 35 26; 36 12; 36 19; 36 37; 36 40; 37 2; 37 5; 37 14; 37 27; 38 12; 38 14; 38 19; 38 27; 38 34; 39 19; 39 30; 39 34; 39 35]
global d_x = [9.0, 3.0, 7.0, 10.0, 9.0, 8.0, 5.0, 5.0, 10.0, 10.0, 7.0, 6.0, 9.0, 2.0, 8.0, 6.0, 9.0, 10.0, 3.0, 5.0, 6.0, 7.0, 3.0, 8.0, 7.0, 3.0, 5.0, 7.0, 2.0, 4.0, 6.0, 10.0, 10.0, 10.0, 3.0, 10.0, 3.0, 9.0, 10.0, 6.0, 8.0, 4.0, 2.0, 3.0, 8.0, 8.0, 6.0, 1.0, 7.0, 9.0, 10.0, 9.0, 7.0, 5.0, 1.0, 5.0, 10.0, 10.0, 7.0, 2.0, 7.0, 5.0, 6.0, 10.0, 7.0, 8.0, 8.0, 1.0, 2.0, 7.0, 3.0, 2.0, 10.0, 10.0, 7.0, 7.0, 5.0, 2.0, 5.0, 4.0, 1.0, 9.0, 6.0, 5.0, 10.0, 1.0, 9.0, 8.0, 7.0, 6.0, 4.0, 2.0, 10.0, 1.0, 8.0, 6.0, 8.0, 2.0, 9.0, 10.0, 2.0, 10.0, 1.0, 7.0, 3.0, 3.0, 10.0, 7.0, 9.0, 7.0, 2.0, 8.0, 3.0, 4.0, 6.0, 2.0, 1.0, 9.0, 8.0, 10.0, 1.0, 9.0, 10.0, 9.0, 8.0, 6.0, 4.0, 2.0, 8.0, 8.0, 2.0, 8.0, 9.0, 10.0, 5.0, 8.0, 1.0, 2.0, 6.0, 5.0, 7.0, 7.0, 6.0]
global b_x = 5
global d_y = [4.0, 7.0, 4.0, 5.0, 10.0, 4.0, 9.0, 4.0, 8.0, 3.0, 3.0, 1.0, 5.0, 6.0, 10.0, 5.0, 2.0, 8.0, 9.0, 7.0, 4.0, 4.0, 2.0, 10.0, 10.0, 6.0, 10.0, 6.0, 1.0, 3.0, 7.0, 6.0, 9.0, 5.0, 2.0, 3.0, 3.0, 5.0, 2.0, 6.0, 9.0, 6.0, 3.0, 6.0, 2.0, 4.0, 5.0, 1.0, 2.0, 1.0, 1.0, 2.0, 4.0, 1.0, 3.0, 1.0, 7.0, 10.0, 5.0, 1.0, 5.0, 9.0, 10.0, 1.0, 8.0, 4.0, 6.0, 9.0, 10.0, 7.0, 3.0, 8.0, 9.0, 8.0, 4.0, 3.0, 10.0, 7.0, 9.0, 3.0, 6.0, 9.0, 4.0, 2.0, 10.0, 2.0, 2.0, 6.0, 7.0, 4.0, 1.0, 2.0, 10.0, 1.0, 2.0, 6.0, 8.0, 3.0, 5.0, 6.0, 5.0, 7.0, 4.0, 6.0, 2.0, 10.0, 8.0, 3.0, 8.0, 6.0, 9.0, 4.0, 3.0, 3.0, 9.0, 9.0, 8.0, 4.0, 6.0, 2.0, 4.0, 8.0, 2.0, 8.0, 2.0, 6.0, 6.0, 3.0, 8.0, 10.0, 3.0, 1.0, 1.0, 4.0, 4.0, 7.0, 5.0, 1.0, 1.0, 9.0, 7.0, 7.0, 8.0]
global b_y = 10
global p = [0.908, 0.828, 0.773, 0.572, 0.098, 0.502, 0.774, 0.584, 0.753, 0.456, 0.85, 0.191, 0.147, 0.704, 0.338, 0.587, 0.481, 0.233, 0.623, 0.829, 0.366, 0.305, 0.639, 0.666, 0.467, 0.316, 0.316, 0.176, 0.132, 0.286, 0.121, 0.919, 0.238, 0.948, 0.719, 0.856, 0.548, 0.021, 0.293, 0.255, 0.705, 0.122, 0.712, 0.009, 0.555, 0.368, 0.972, 0.472, 0.527, 0.936, 0.305, 0.864, 0.925, 0.816, 0.936, 0.335, 0.36, 0.737, 0.4, 0.517, 0.872, 0.46, 0.101, 0.458, 0.783, 0.728, 0.37, 0.386, 0.591, 0.813, 0.085, 0.125, 0.306, 0.722, 0.605, 0.001, 0.325, 0.947, 0.523, 0.398, 0.951, 0.967, 0.466, 0.339, 0.201, 0.819, 0.99, 0.519, 0.186, 0.427, 0.691, 0.875, 0.811, 0.482, 0.516, 0.027, 0.164, 0.159, 0.566, 0.715, 0.716, 0.803, 0.156, 0.437, 0.597, 0.035, 0.987, 0.552, 0.08, 0.729, 0.924, 0.464, 0.444, 0.127, 0.768, 0.507, 0.609, 0.1, 0.423, 0.942, 0.507, 0.959, 0.891, 0.919, 0.102, 0.769, 0.171, 0.897, 0.343, 0.919, 0.278, 0.642, 0.84, 0.454, 0.488, 0.109, 0.779, 0.644, 0.968, 0.672, 0.654, 0.66, 0.863]
global q = [0.977, 0.952, 0.874, 0.811, 0.906, 0.695, 0.968, 0.956, 0.992, 0.912, 0.988, 0.248, 0.379, 0.895, 0.725, 0.873, 0.876, 0.367, 0.82, 0.925, 0.693, 0.903, 0.776, 0.862, 0.83, 0.684, 0.848, 0.934, 0.734, 0.814, 0.353, 0.993, 0.945, 0.953, 0.809, 0.909, 0.8, 0.513, 0.347, 0.755, 0.873, 0.884, 0.949, 0.255, 0.845, 0.986, 0.999, 0.856, 0.775, 0.993, 0.804, 0.879, 0.971, 0.9, 0.994, 0.78, 0.755, 0.785, 0.745, 0.916, 0.99, 0.767, 0.13, 0.655, 0.815, 0.769, 0.571, 0.933, 0.931, 0.977, 0.222, 0.978, 0.826, 0.797, 0.951, 0.772, 0.694, 0.967, 0.863, 0.544, 0.952, 0.972, 0.737, 0.608, 0.999, 0.959, 0.996, 0.758, 0.477, 0.629, 0.761, 0.945, 0.919, 0.616, 0.573, 0.24, 0.368, 0.275, 0.725, 0.932, 0.996, 0.909, 0.818, 0.822, 0.916, 0.623, 0.994, 0.94, 0.637, 0.865, 0.937, 0.934, 0.693, 0.297, 0.779, 0.813, 0.898, 0.992, 0.482, 0.979, 0.555, 0.967, 0.938, 0.981, 0.362, 0.858, 0.678, 0.959, 0.691, 0.99, 0.795, 0.687, 0.946, 0.77, 0.755, 0.55, 0.814, 0.756, 0.997, 0.925, 0.779, 0.862, 0.959]
global origin = 1
global destination = 40