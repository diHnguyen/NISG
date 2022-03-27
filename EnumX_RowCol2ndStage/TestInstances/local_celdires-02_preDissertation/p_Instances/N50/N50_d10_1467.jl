global arcs = [1 7; 1 10; 1 16; 1 18; 1 27; 1 37; 2 5; 2 14; 2 28; 2 33; 2 45; 2 47; 2 50; 3 14; 3 20; 3 33; 3 37; 3 38; 4 5; 4 12; 4 33; 4 39; 4 40; 4 41; 4 47; 4 48; 5 7; 5 15; 5 22; 5 28; 5 46; 5 49; 6 18; 6 28; 6 29; 6 33; 6 37; 7 9; 7 17; 7 25; 7 38; 7 42; 7 45; 8 14; 8 17; 8 26; 9 10; 9 33; 9 47; 9 49; 10 27; 11 13; 11 35; 11 50; 12 19; 12 49; 13 3; 13 10; 13 12; 13 14; 13 24; 13 32; 13 40; 13 47; 13 50; 14 31; 14 47; 15 6; 15 9; 15 17; 15 19; 16 4; 16 28; 16 29; 16 34; 16 42; 17 11; 17 26; 17 28; 18 20; 18 39; 18 42; 19 8; 19 17; 19 20; 19 22; 19 31; 19 33; 19 44; 19 47; 20 12; 20 16; 20 19; 20 27; 20 45; 20 46; 20 50; 21 2; 21 4; 21 11; 21 12; 21 14; 21 26; 21 38; 22 2; 22 6; 22 8; 22 16; 22 32; 23 11; 23 24; 24 7; 24 11; 24 16; 24 22; 24 25; 24 31; 24 49; 25 13; 25 16; 25 29; 26 24; 26 27; 26 30; 26 35; 26 40; 26 41; 27 15; 27 34; 27 46; 28 6; 28 7; 28 9; 28 21; 28 39; 29 21; 29 47; 30 2; 30 6; 30 8; 30 11; 30 12; 30 31; 30 41; 30 43; 30 47; 31 26; 31 29; 31 35; 31 44; 31 46; 32 30; 32 34; 32 35; 32 36; 32 41; 32 42; 33 4; 33 6; 33 28; 33 40; 33 48; 34 4; 34 19; 34 23; 34 24; 34 26; 34 49; 35 26; 35 48; 36 14; 36 37; 36 48; 37 8; 37 10; 37 25; 37 27; 37 30; 37 32; 37 38; 38 5; 38 8; 38 17; 38 28; 38 35; 38 39; 38 43; 38 50; 39 19; 39 29; 40 9; 40 13; 40 36; 40 39; 41 4; 41 16; 41 27; 41 33; 42 12; 42 28; 42 29; 42 39; 42 44; 43 24; 43 33; 43 46; 44 11; 45 3; 45 21; 45 22; 45 30; 46 5; 46 32; 47 12; 47 17; 47 22; 47 27; 47 39; 48 3; 48 5; 48 7; 48 8; 48 12; 48 16; 48 25; 48 26; 48 29; 48 32; 48 40; 48 42; 49 16; 49 19; 49 35; 49 36; 49 37; 49 43; 49 45]
global d_x = [8.0, 9.0, 5.0, 9.0, 4.0, 7.0, 9.0, 10.0, 6.0, 7.0, 10.0, 1.0, 5.0, 5.0, 5.0, 4.0, 4.0, 1.0, 6.0, 9.0, 7.0, 3.0, 7.0, 3.0, 1.0, 10.0, 1.0, 5.0, 2.0, 6.0, 6.0, 9.0, 1.0, 1.0, 6.0, 1.0, 1.0, 8.0, 2.0, 5.0, 5.0, 5.0, 7.0, 8.0, 9.0, 10.0, 1.0, 5.0, 10.0, 1.0, 6.0, 4.0, 3.0, 5.0, 8.0, 2.0, 2.0, 9.0, 6.0, 4.0, 2.0, 2.0, 5.0, 4.0, 2.0, 10.0, 9.0, 10.0, 6.0, 1.0, 9.0, 8.0, 2.0, 5.0, 2.0, 1.0, 6.0, 8.0, 7.0, 7.0, 2.0, 1.0, 6.0, 8.0, 4.0, 1.0, 8.0, 1.0, 3.0, 7.0, 1.0, 2.0, 1.0, 1.0, 6.0, 1.0, 7.0, 8.0, 9.0, 8.0, 9.0, 2.0, 5.0, 1.0, 8.0, 4.0, 4.0, 6.0, 8.0, 5.0, 7.0, 9.0, 8.0, 9.0, 9.0, 9.0, 5.0, 3.0, 5.0, 10.0, 2.0, 10.0, 9.0, 7.0, 2.0, 5.0, 10.0, 9.0, 2.0, 8.0, 9.0, 3.0, 5.0, 1.0, 7.0, 1.0, 5.0, 5.0, 10.0, 9.0, 2.0, 5.0, 3.0, 10.0, 7.0, 6.0, 5.0, 1.0, 4.0, 4.0, 9.0, 6.0, 9.0, 2.0, 7.0, 6.0, 3.0, 5.0, 7.0, 9.0, 8.0, 4.0, 8.0, 2.0, 6.0, 7.0, 6.0, 8.0, 4.0, 6.0, 4.0, 10.0, 3.0, 7.0, 2.0, 8.0, 7.0, 9.0, 10.0, 5.0, 10.0, 4.0, 1.0, 7.0, 4.0, 10.0, 9.0, 4.0, 6.0, 9.0, 9.0, 3.0, 1.0, 2.0, 7.0, 1.0, 8.0, 1.0, 8.0, 2.0, 6.0, 9.0, 8.0, 10.0, 8.0, 10.0, 10.0, 10.0, 5.0, 10.0, 8.0, 9.0, 1.0, 3.0, 4.0, 2.0, 8.0, 5.0, 1.0, 7.0, 5.0, 2.0, 8.0, 10.0, 8.0, 3.0, 8.0, 6.0, 10.0, 7.0, 3.0, 7.0, 7.0, 7.0, 8.0, 8.0, 3.0]
global b_x = 5
global d_y = [7.0, 4.0, 8.0, 5.0, 2.0, 7.0, 7.0, 4.0, 5.0, 2.0, 8.0, 9.0, 5.0, 9.0, 4.0, 7.0, 4.0, 8.0, 1.0, 2.0, 3.0, 4.0, 8.0, 3.0, 6.0, 5.0, 8.0, 4.0, 7.0, 1.0, 9.0, 6.0, 9.0, 1.0, 3.0, 10.0, 8.0, 6.0, 2.0, 7.0, 7.0, 1.0, 10.0, 8.0, 3.0, 7.0, 5.0, 8.0, 3.0, 5.0, 10.0, 8.0, 8.0, 2.0, 10.0, 5.0, 10.0, 4.0, 10.0, 10.0, 10.0, 10.0, 9.0, 6.0, 7.0, 9.0, 7.0, 6.0, 5.0, 10.0, 2.0, 1.0, 9.0, 2.0, 3.0, 9.0, 9.0, 1.0, 3.0, 3.0, 7.0, 2.0, 9.0, 2.0, 8.0, 1.0, 4.0, 10.0, 6.0, 7.0, 5.0, 9.0, 6.0, 9.0, 6.0, 4.0, 2.0, 3.0, 5.0, 5.0, 1.0, 6.0, 2.0, 5.0, 1.0, 5.0, 4.0, 1.0, 8.0, 6.0, 6.0, 2.0, 7.0, 4.0, 10.0, 6.0, 5.0, 2.0, 7.0, 9.0, 10.0, 1.0, 5.0, 4.0, 5.0, 6.0, 10.0, 9.0, 2.0, 5.0, 1.0, 9.0, 5.0, 9.0, 10.0, 2.0, 5.0, 9.0, 7.0, 6.0, 10.0, 7.0, 2.0, 5.0, 7.0, 2.0, 10.0, 6.0, 6.0, 6.0, 4.0, 7.0, 10.0, 3.0, 3.0, 1.0, 1.0, 9.0, 2.0, 10.0, 6.0, 5.0, 6.0, 3.0, 8.0, 8.0, 4.0, 7.0, 1.0, 2.0, 6.0, 1.0, 1.0, 8.0, 3.0, 5.0, 5.0, 2.0, 4.0, 3.0, 6.0, 3.0, 7.0, 6.0, 7.0, 1.0, 9.0, 5.0, 7.0, 10.0, 10.0, 10.0, 4.0, 3.0, 2.0, 8.0, 9.0, 6.0, 9.0, 3.0, 10.0, 1.0, 6.0, 4.0, 8.0, 1.0, 7.0, 10.0, 9.0, 3.0, 2.0, 3.0, 7.0, 2.0, 7.0, 5.0, 1.0, 10.0, 10.0, 4.0, 9.0, 2.0, 5.0, 8.0, 4.0, 2.0, 4.0, 9.0, 10.0, 4.0, 4.0, 1.0, 9.0, 2.0, 2.0, 2.0, 2.0]
global b_y = 10
global p = [0.617, 0.674, 0.35, 0.772, 0.95, 0.873, 0.776, 0.577, 0.112, 0.654, 0.552, 0.817, 0.542, 0.117, 0.282, 0.025, 0.814, 0.812, 0.723, 0.148, 0.919, 0.305, 0.99, 0.861, 0.49, 0.155, 0.573, 0.806, 0.305, 0.204, 0.014, 0.655, 0.525, 0.01, 0.604, 0.543, 0.667, 0.559, 0.471, 0.449, 0.729, 0.435, 0.951, 0.724, 0.163, 0.044, 0.796, 0.906, 0.23, 0.297, 0.717, 0.858, 0.025, 0.205, 0.271, 0.726, 0.639, 0.558, 0.829, 0.224, 0.052, 0.26, 0.067, 0.591, 0.593, 0.866, 0.577, 0.338, 0.344, 0.087, 0.143, 0.842, 0.414, 0.327, 0.388, 0.614, 0.292, 0.811, 0.992, 0.684, 0.372, 0.766, 0.738, 0.133, 0.266, 0.083, 0.857, 0.598, 0.046, 0.768, 0.728, 0.561, 0.015, 0.426, 0.47, 0.965, 0.11, 0.473, 0.563, 0.366, 0.64, 0.078, 0.716, 0.959, 0.664, 0.939, 0.421, 0.772, 0.498, 0.963, 0.231, 0.82, 0.661, 0.995, 0.585, 0.682, 0.107, 0.19, 0.737, 0.908, 0.056, 0.054, 0.758, 0.304, 0.935, 0.27, 0.063, 0.754, 0.026, 0.042, 0.754, 0.579, 0.906, 0.626, 0.39, 0.531, 0.188, 0.366, 0.991, 0.37, 0.139, 0.552, 0.704, 0.853, 0.362, 0.146, 0.304, 0.244, 0.616, 0.959, 0.623, 0.39, 0.905, 0.381, 0.463, 0.183, 0.181, 0.466, 0.506, 0.086, 0.238, 0.172, 0.14, 0.088, 0.522, 0.392, 0.017, 0.056, 0.046, 0.257, 0.424, 0.618, 0.814, 0.746, 0.445, 0.616, 0.99, 0.058, 0.864, 0.481, 0.307, 0.64, 0.625, 0.748, 0.266, 0.859, 0.042, 0.46, 0.051, 0.583, 0.059, 0.826, 0.909, 0.005, 0.266, 0.706, 0.095, 0.753, 0.463, 0.15, 0.215, 0.678, 0.189, 0.755, 0.247, 0.096, 0.685, 0.283, 0.055, 0.674, 0.106, 0.122, 0.518, 0.246, 0.451, 0.508, 0.277, 0.82, 0.152, 0.765, 0.42, 0.428, 0.513, 0.113, 0.473, 0.16, 0.272, 0.707, 0.956, 0.556, 0.621, 0.394, 0.917, 0.598, 0.815, 0.41, 0.474]
global q = [0.879, 0.978, 0.544, 0.959, 0.981, 0.981, 0.87, 0.607, 0.802, 0.674, 0.63, 0.974, 0.869, 0.603, 0.286, 0.468, 0.873, 0.946, 0.832, 0.842, 0.926, 0.318, 0.994, 0.862, 0.662, 0.85, 0.94, 0.818, 0.924, 0.793, 0.93, 0.656, 0.591, 0.662, 0.751, 0.747, 0.737, 0.649, 0.73, 0.449, 0.908, 0.53, 0.978, 0.737, 0.565, 0.247, 0.99, 0.99, 0.883, 0.519, 0.811, 0.985, 0.439, 0.559, 0.992, 0.886, 0.731, 0.746, 0.848, 0.36, 0.112, 0.737, 0.864, 0.857, 0.941, 0.866, 0.591, 0.52, 0.443, 0.489, 0.997, 0.981, 0.593, 0.524, 0.81, 0.734, 0.735, 0.987, 0.998, 0.979, 0.921, 0.878, 0.956, 0.644, 0.28, 0.536, 0.86, 0.885, 0.174, 0.821, 0.902, 0.639, 0.26, 0.523, 0.507, 0.971, 0.87, 0.532, 0.742, 0.691, 0.698, 0.487, 0.903, 0.99, 0.898, 0.97, 0.525, 0.875, 0.642, 0.985, 0.296, 0.975, 0.751, 0.996, 0.714, 0.876, 0.199, 0.745, 0.75, 0.981, 0.584, 0.512, 0.791, 0.907, 0.943, 0.518, 0.465, 0.843, 0.512, 0.859, 0.995, 0.884, 0.909, 0.658, 0.653, 0.883, 0.509, 0.783, 0.994, 0.409, 0.3, 0.677, 0.956, 0.895, 0.642, 0.736, 0.811, 0.634, 0.702, 0.975, 0.997, 0.67, 0.966, 0.921, 0.676, 0.469, 0.833, 0.954, 0.674, 0.985, 0.662, 0.379, 0.605, 0.982, 0.887, 0.717, 0.744, 0.13, 0.24, 0.591, 0.606, 0.797, 0.978, 0.835, 0.846, 0.636, 0.99, 0.206, 0.931, 0.831, 0.73, 0.64, 0.664, 0.803, 0.852, 0.968, 0.75, 0.651, 0.99, 0.711, 0.352, 0.852, 0.949, 0.148, 0.396, 0.841, 0.622, 0.994, 0.589, 0.19, 0.399, 0.786, 0.682, 0.775, 0.883, 0.591, 0.811, 0.787, 0.693, 0.799, 0.566, 0.735, 0.943, 0.79, 0.904, 0.945, 0.435, 0.827, 0.43, 0.979, 0.605, 0.87, 0.954, 0.218, 0.857, 0.202, 0.399, 0.853, 0.957, 0.776, 0.745, 0.818, 0.949, 0.994, 0.827, 0.609, 0.791]
global origin = 1
global destination = 50