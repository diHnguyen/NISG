global arcs = [1 13; 1 20; 1 21; 1 24; 1 25; 1 33; 1 39; 1 40; 1 49; 2 14; 2 36; 2 44; 2 48; 3 7; 3 19; 3 22; 3 31; 3 33; 3 36; 3 37; 3 41; 4 17; 4 30; 4 36; 4 38; 4 39; 4 41; 5 7; 5 9; 5 18; 5 42; 5 47; 5 48; 5 50; 6 10; 6 39; 6 46; 7 3; 7 23; 7 26; 7 32; 7 45; 7 49; 8 3; 8 10; 8 20; 8 22; 8 34; 8 40; 8 42; 9 4; 9 7; 9 23; 9 30; 9 36; 10 9; 10 11; 10 24; 10 45; 10 49; 11 6; 11 7; 11 21; 11 29; 11 31; 11 40; 12 5; 12 31; 12 33; 12 37; 12 41; 13 14; 13 16; 13 19; 13 24; 13 25; 13 43; 13 46; 14 36; 14 37; 14 50; 15 9; 15 10; 15 29; 15 35; 15 36; 15 37; 15 45; 16 4; 16 29; 16 30; 16 47; 16 49; 17 2; 17 13; 17 30; 18 19; 18 31; 18 32; 18 35; 19 3; 19 4; 19 21; 19 43; 19 47; 20 15; 20 26; 20 33; 20 40; 20 41; 20 50; 21 9; 21 12; 21 40; 21 46; 22 5; 22 6; 22 23; 22 29; 22 41; 22 46; 23 4; 23 7; 23 18; 23 20; 23 21; 23 34; 23 40; 23 43; 23 48; 24 7; 24 8; 24 11; 24 15; 24 17; 24 18; 24 30; 24 36; 25 22; 26 14; 26 30; 26 43; 27 10; 27 23; 28 49; 29 3; 29 7; 29 9; 29 20; 29 41; 29 42; 29 45; 30 9; 30 25; 30 29; 30 37; 30 45; 30 50; 31 16; 31 33; 31 39; 31 43; 31 50; 32 16; 32 26; 32 27; 32 33; 32 45; 33 16; 33 19; 33 45; 34 33; 34 40; 34 42; 35 7; 35 10; 35 32; 35 39; 35 41; 35 47; 36 3; 36 40; 37 32; 37 38; 37 48; 38 4; 38 15; 38 43; 39 2; 39 32; 39 35; 39 40; 39 42; 40 12; 40 14; 40 21; 40 42; 40 44; 40 48; 41 5; 41 7; 41 43; 41 44; 42 8; 42 46; 43 4; 43 18; 43 20; 43 28; 43 30; 43 31; 43 35; 43 41; 43 46; 44 10; 44 17; 44 23; 44 27; 44 40; 45 3; 45 5; 45 25; 45 28; 45 30; 45 34; 45 35; 46 13; 46 14; 46 21; 46 25; 46 35; 46 44; 47 4; 47 5; 47 7; 47 15; 47 16; 47 30; 47 36; 47 40; 47 45; 48 8; 48 19; 48 24; 48 26; 48 44; 48 45; 48 50; 49 21; 49 29; 49 36]
global d_x = [3.0, 7.0, 3.0, 3.0, 9.0, 6.0, 7.0, 1.0, 7.0, 1.0, 9.0, 9.0, 9.0, 10.0, 3.0, 3.0, 9.0, 7.0, 4.0, 3.0, 4.0, 6.0, 7.0, 1.0, 5.0, 8.0, 3.0, 2.0, 7.0, 3.0, 7.0, 9.0, 5.0, 6.0, 10.0, 4.0, 3.0, 5.0, 4.0, 7.0, 3.0, 7.0, 2.0, 10.0, 9.0, 10.0, 3.0, 1.0, 1.0, 4.0, 7.0, 10.0, 3.0, 3.0, 4.0, 4.0, 5.0, 4.0, 1.0, 10.0, 6.0, 10.0, 2.0, 7.0, 8.0, 7.0, 2.0, 8.0, 7.0, 9.0, 8.0, 9.0, 3.0, 4.0, 6.0, 1.0, 1.0, 5.0, 7.0, 3.0, 5.0, 8.0, 5.0, 2.0, 8.0, 2.0, 2.0, 3.0, 8.0, 10.0, 5.0, 1.0, 1.0, 8.0, 9.0, 8.0, 4.0, 8.0, 4.0, 8.0, 3.0, 8.0, 9.0, 5.0, 7.0, 5.0, 9.0, 5.0, 10.0, 4.0, 8.0, 10.0, 1.0, 8.0, 8.0, 8.0, 10.0, 5.0, 2.0, 7.0, 5.0, 6.0, 3.0, 4.0, 9.0, 10.0, 9.0, 8.0, 8.0, 1.0, 6.0, 7.0, 6.0, 8.0, 1.0, 1.0, 1.0, 10.0, 4.0, 4.0, 9.0, 6.0, 4.0, 3.0, 7.0, 6.0, 5.0, 7.0, 3.0, 6.0, 5.0, 2.0, 5.0, 9.0, 7.0, 9.0, 3.0, 2.0, 9.0, 4.0, 8.0, 2.0, 5.0, 1.0, 1.0, 5.0, 2.0, 4.0, 5.0, 2.0, 1.0, 1.0, 9.0, 5.0, 3.0, 8.0, 6.0, 9.0, 5.0, 3.0, 9.0, 2.0, 7.0, 4.0, 6.0, 3.0, 4.0, 5.0, 4.0, 2.0, 5.0, 7.0, 10.0, 8.0, 5.0, 5.0, 7.0, 8.0, 10.0, 9.0, 2.0, 4.0, 6.0, 8.0, 4.0, 10.0, 4.0, 2.0, 7.0, 1.0, 6.0, 9.0, 1.0, 9.0, 4.0, 1.0, 10.0, 5.0, 6.0, 4.0, 3.0, 3.0, 3.0, 6.0, 4.0, 7.0, 2.0, 1.0, 2.0, 8.0, 5.0, 5.0, 2.0, 4.0, 7.0, 1.0, 4.0, 10.0, 3.0, 2.0, 8.0, 8.0, 2.0, 1.0, 1.0, 2.0, 4.0, 6.0, 8.0, 5.0, 10.0]
global b_x = 5
global d_y = [10.0, 9.0, 7.0, 6.0, 9.0, 9.0, 9.0, 6.0, 6.0, 1.0, 6.0, 1.0, 9.0, 10.0, 1.0, 7.0, 6.0, 7.0, 2.0, 2.0, 3.0, 1.0, 4.0, 4.0, 8.0, 6.0, 9.0, 7.0, 1.0, 2.0, 4.0, 2.0, 7.0, 3.0, 4.0, 5.0, 1.0, 9.0, 7.0, 3.0, 2.0, 2.0, 7.0, 4.0, 5.0, 1.0, 6.0, 1.0, 10.0, 1.0, 3.0, 7.0, 4.0, 9.0, 5.0, 1.0, 2.0, 7.0, 8.0, 4.0, 8.0, 2.0, 4.0, 3.0, 7.0, 10.0, 4.0, 7.0, 10.0, 3.0, 7.0, 3.0, 1.0, 5.0, 6.0, 2.0, 7.0, 8.0, 10.0, 7.0, 3.0, 6.0, 3.0, 5.0, 4.0, 8.0, 3.0, 7.0, 9.0, 5.0, 9.0, 3.0, 5.0, 5.0, 10.0, 10.0, 10.0, 6.0, 5.0, 7.0, 4.0, 8.0, 3.0, 2.0, 7.0, 7.0, 2.0, 5.0, 10.0, 6.0, 8.0, 5.0, 7.0, 9.0, 3.0, 8.0, 2.0, 1.0, 4.0, 4.0, 10.0, 7.0, 7.0, 3.0, 9.0, 7.0, 4.0, 9.0, 6.0, 3.0, 4.0, 10.0, 9.0, 4.0, 9.0, 1.0, 4.0, 1.0, 10.0, 3.0, 1.0, 4.0, 7.0, 2.0, 3.0, 2.0, 9.0, 2.0, 2.0, 9.0, 6.0, 6.0, 2.0, 6.0, 2.0, 8.0, 2.0, 2.0, 5.0, 9.0, 6.0, 9.0, 8.0, 2.0, 2.0, 6.0, 10.0, 5.0, 10.0, 10.0, 2.0, 10.0, 5.0, 3.0, 1.0, 10.0, 6.0, 6.0, 10.0, 4.0, 8.0, 1.0, 10.0, 1.0, 4.0, 1.0, 7.0, 10.0, 6.0, 2.0, 7.0, 8.0, 5.0, 8.0, 10.0, 5.0, 4.0, 3.0, 3.0, 9.0, 8.0, 2.0, 5.0, 2.0, 1.0, 8.0, 4.0, 7.0, 1.0, 1.0, 10.0, 4.0, 10.0, 6.0, 8.0, 10.0, 3.0, 6.0, 10.0, 5.0, 6.0, 1.0, 4.0, 10.0, 8.0, 4.0, 3.0, 9.0, 9.0, 3.0, 1.0, 1.0, 7.0, 4.0, 1.0, 6.0, 3.0, 6.0, 2.0, 10.0, 9.0, 5.0, 2.0, 1.0, 10.0, 2.0, 4.0, 2.0, 2.0, 4.0, 10.0]
global b_y = 10
global p = [0.825, 0.164, 0.17, 0.182, 0.312, 0.854, 0.224, 0.661, 0.635, 0.607, 0.09, 0.249, 0.283, 0.119, 0.652, 0.719, 0.855, 0.941, 0.082, 0.818, 0.445, 0.853, 0.677, 0.16, 0.864, 0.208, 0.092, 0.878, 0.504, 0.372, 0.863, 0.214, 0.648, 0.201, 0.147, 0.574, 0.771, 0.723, 0.079, 0.01, 0.571, 0.364, 0.573, 0.681, 0.635, 0.688, 0.368, 0.261, 0.431, 0.427, 0.268, 0.597, 0.887, 0.953, 0.303, 0.963, 0.655, 0.513, 0.981, 0.634, 0.869, 0.44, 0.966, 0.272, 0.207, 0.947, 0.93, 0.819, 0.781, 0.532, 0.613, 0.781, 0.486, 0.273, 0.467, 0.454, 0.562, 0.779, 0.317, 0.049, 0.803, 0.397, 0.248, 0.543, 0.026, 0.769, 0.549, 0.92, 0.926, 0.733, 0.102, 0.182, 0.935, 0.623, 0.203, 0.421, 0.421, 0.338, 0.446, 0.448, 0.021, 0.394, 0.586, 0.401, 0.546, 0.256, 0.949, 0.748, 0.552, 0.575, 0.593, 0.333, 0.432, 0.762, 0.784, 0.235, 0.585, 0.052, 0.159, 0.925, 0.034, 0.628, 0.011, 0.887, 0.736, 0.862, 0.243, 0.441, 0.123, 0.601, 0.555, 0.285, 0.051, 0.281, 0.502, 0.608, 0.775, 0.553, 0.94, 0.227, 0.626, 0.608, 0.062, 0.07, 0.843, 0.005, 0.984, 0.797, 0.747, 0.861, 0.409, 0.132, 0.3, 0.874, 0.01, 0.765, 0.727, 0.151, 0.56, 0.414, 0.218, 0.854, 0.826, 0.969, 0.732, 0.495, 0.183, 0.678, 0.903, 0.436, 0.391, 0.461, 0.277, 0.01, 0.924, 0.768, 0.717, 0.983, 0.985, 0.691, 0.111, 0.011, 0.857, 0.328, 0.528, 0.663, 0.119, 0.437, 0.098, 0.058, 0.275, 0.826, 0.882, 0.414, 0.698, 0.263, 0.65, 0.23, 0.798, 0.708, 0.488, 0.174, 0.328, 0.757, 0.276, 0.375, 0.976, 0.623, 0.211, 0.306, 0.59, 0.665, 0.203, 0.74, 0.73, 0.176, 0.996, 0.392, 0.027, 0.701, 0.263, 0.118, 0.366, 0.802, 0.799, 0.438, 0.916, 0.514, 0.685, 0.336, 0.461, 0.932, 0.066, 0.049, 0.363, 0.806, 0.929, 0.967, 0.717, 0.498, 0.596, 0.894, 0.541, 0.437, 0.314, 0.655, 0.368, 0.121, 0.012, 0.12, 0.649]
global q = [0.88, 0.871, 0.215, 0.457, 0.365, 0.972, 0.91, 0.803, 0.886, 0.746, 0.658, 0.998, 0.603, 0.563, 0.772, 0.759, 0.969, 0.945, 0.089, 0.995, 0.477, 0.918, 0.969, 0.168, 0.89, 0.409, 0.839, 0.908, 0.61, 0.977, 0.938, 0.441, 0.863, 0.481, 0.149, 0.781, 0.861, 0.804, 0.862, 0.907, 0.591, 0.918, 0.594, 0.758, 0.727, 0.841, 0.515, 0.294, 0.933, 0.638, 0.316, 0.689, 0.993, 0.956, 0.4, 0.986, 0.963, 0.526, 0.983, 0.708, 0.932, 0.533, 0.986, 0.386, 0.342, 0.989, 0.954, 0.845, 0.841, 0.755, 0.753, 0.969, 0.996, 0.611, 0.619, 0.507, 0.973, 0.816, 0.925, 0.143, 0.965, 0.401, 0.697, 0.641, 0.748, 0.839, 0.662, 0.958, 0.955, 0.778, 0.453, 0.278, 0.95, 0.813, 0.898, 0.791, 0.986, 0.652, 0.655, 0.472, 0.112, 0.915, 0.775, 0.66, 0.577, 0.67, 0.979, 0.82, 0.875, 0.645, 0.631, 0.392, 0.522, 0.877, 0.88, 0.505, 0.637, 0.938, 0.358, 0.932, 0.738, 0.973, 0.274, 0.997, 0.823, 0.926, 0.797, 0.598, 0.24, 0.828, 0.588, 0.647, 0.632, 0.308, 0.773, 0.906, 0.902, 0.847, 0.982, 0.479, 0.688, 0.896, 0.46, 0.48, 0.972, 0.06, 0.997, 0.999, 0.774, 0.979, 0.503, 0.195, 0.643, 0.935, 0.97, 0.79, 0.98, 0.165, 0.91, 0.81, 0.4, 0.955, 0.83, 0.993, 0.943, 0.776, 0.693, 0.966, 0.919, 0.494, 0.614, 0.494, 0.628, 0.433, 0.968, 0.797, 0.777, 0.988, 0.993, 0.704, 0.432, 0.67, 0.955, 0.846, 0.535, 0.942, 0.691, 0.813, 0.778, 0.489, 0.789, 0.911, 0.899, 0.538, 0.826, 0.5, 0.8, 0.654, 0.823, 0.931, 0.898, 0.367, 0.447, 0.986, 0.767, 0.986, 0.98, 0.747, 0.612, 0.846, 0.682, 0.7, 0.731, 0.909, 0.837, 0.461, 0.998, 0.439, 0.562, 0.762, 0.345, 0.858, 0.515, 0.848, 0.814, 0.936, 0.92, 0.987, 0.962, 0.821, 0.815, 0.999, 0.384, 0.204, 0.689, 0.814, 0.948, 0.986, 0.784, 0.657, 0.645, 0.923, 0.935, 0.509, 0.917, 0.673, 0.941, 0.201, 0.311, 0.203, 0.891]
global origin = 1
global destination = 50