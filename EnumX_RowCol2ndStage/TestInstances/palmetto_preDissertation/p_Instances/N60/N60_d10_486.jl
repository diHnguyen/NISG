global arcs = [1 10; 1 16; 1 24; 1 39; 1 55; 1 56; 1 57; 2 8; 2 14; 2 17; 2 42; 3 9; 3 21; 3 59; 4 8; 4 17; 4 23; 4 24; 4 26; 4 36; 4 47; 4 49; 5 40; 5 42; 5 45; 6 3; 6 13; 6 28; 6 34; 6 57; 6 59; 6 60; 7 13; 7 33; 7 38; 7 56; 8 44; 8 46; 8 56; 8 60; 9 2; 9 8; 9 10; 9 16; 9 43; 9 48; 10 12; 10 16; 10 18; 10 29; 10 41; 10 50; 10 59; 11 5; 11 8; 11 9; 11 15; 11 16; 11 18; 11 24; 11 48; 11 60; 12 3; 12 9; 12 14; 12 19; 12 29; 12 53; 13 6; 13 26; 13 31; 13 36; 13 37; 13 38; 13 46; 14 3; 14 7; 14 20; 14 27; 14 34; 14 38; 15 13; 15 18; 15 36; 15 43; 15 48; 16 4; 16 8; 16 15; 16 21; 16 41; 16 44; 16 46; 16 53; 16 57; 17 11; 17 31; 17 41; 17 58; 18 4; 18 12; 18 13; 18 21; 18 27; 18 36; 18 46; 18 50; 18 51; 19 9; 19 10; 19 25; 19 27; 19 42; 19 57; 20 2; 20 10; 20 14; 20 28; 20 41; 20 44; 20 46; 20 48; 21 9; 21 10; 21 20; 21 29; 21 37; 22 30; 22 32; 22 59; 23 5; 23 15; 23 17; 23 45; 23 60; 24 11; 24 20; 24 39; 24 58; 25 19; 25 33; 25 46; 25 50; 26 5; 26 24; 26 39; 26 49; 26 52; 27 13; 27 24; 27 38; 27 49; 27 58; 28 4; 28 6; 28 19; 28 35; 28 38; 28 42; 28 46; 28 48; 28 54; 28 57; 29 2; 29 11; 29 21; 29 31; 29 38; 29 46; 30 4; 30 6; 30 21; 30 22; 30 42; 30 49; 31 12; 31 13; 31 22; 31 39; 31 51; 32 13; 32 35; 33 10; 33 16; 33 31; 33 52; 34 5; 34 11; 34 14; 34 24; 35 13; 35 17; 35 24; 35 38; 35 45; 35 58; 36 21; 36 42; 36 47; 36 50; 36 60; 37 2; 37 7; 37 9; 37 16; 37 31; 37 39; 38 2; 38 36; 38 42; 38 51; 38 55; 38 56; 38 59; 39 13; 39 19; 39 26; 39 52; 39 60; 40 10; 40 12; 40 18; 40 20; 40 26; 40 31; 40 42; 40 60; 41 13; 41 19; 41 20; 41 21; 41 22; 41 25; 41 33; 41 38; 41 52; 42 8; 42 10; 42 17; 42 18; 42 30; 42 47; 43 3; 43 19; 43 31; 43 57; 44 8; 44 11; 44 33; 45 11; 45 17; 45 22; 45 23; 45 34; 45 46; 45 48; 45 53; 46 10; 46 17; 46 28; 46 36; 46 38; 47 6; 47 12; 47 13; 47 28; 47 31; 47 41; 47 42; 47 44; 47 48; 47 52; 47 55; 47 58; 48 17; 48 20; 48 23; 48 30; 48 34; 48 35; 48 37; 48 47; 48 49; 49 8; 49 44; 49 47; 50 9; 50 10; 50 14; 50 16; 50 22; 50 51; 50 56; 50 57; 51 4; 51 26; 51 33; 51 34; 51 35; 51 50; 52 2; 52 4; 52 28; 52 31; 52 49; 53 21; 54 6; 54 28; 54 38; 55 17; 55 57; 55 60; 56 6; 56 22; 56 23; 56 25; 56 28; 56 29; 56 45; 57 35; 57 48; 57 52; 57 54; 57 58; 58 9; 58 11; 58 16; 58 18; 58 19; 58 28; 58 29; 58 30; 58 31; 58 42; 58 55; 59 11; 59 24; 59 26; 59 34; 59 42; 59 53; 59 56]
global d_x = [5.0, 7.0, 4.0, 10.0, 8.0, 8.0, 6.0, 9.0, 9.0, 5.0, 10.0, 8.0, 2.0, 8.0, 8.0, 8.0, 6.0, 5.0, 10.0, 5.0, 10.0, 5.0, 1.0, 6.0, 3.0, 10.0, 3.0, 9.0, 5.0, 8.0, 8.0, 3.0, 1.0, 10.0, 4.0, 8.0, 4.0, 10.0, 6.0, 10.0, 6.0, 6.0, 8.0, 6.0, 4.0, 8.0, 3.0, 2.0, 5.0, 2.0, 8.0, 9.0, 2.0, 10.0, 2.0, 5.0, 10.0, 7.0, 3.0, 8.0, 4.0, 1.0, 3.0, 9.0, 7.0, 9.0, 10.0, 10.0, 6.0, 5.0, 1.0, 5.0, 8.0, 3.0, 4.0, 7.0, 10.0, 1.0, 1.0, 1.0, 9.0, 7.0, 9.0, 10.0, 2.0, 3.0, 10.0, 5.0, 2.0, 7.0, 9.0, 4.0, 9.0, 10.0, 1.0, 2.0, 8.0, 1.0, 9.0, 8.0, 9.0, 8.0, 10.0, 7.0, 8.0, 9.0, 4.0, 6.0, 9.0, 9.0, 10.0, 2.0, 7.0, 6.0, 5.0, 3.0, 7.0, 6.0, 6.0, 3.0, 7.0, 7.0, 3.0, 3.0, 4.0, 1.0, 3.0, 7.0, 10.0, 10.0, 3.0, 4.0, 4.0, 3.0, 1.0, 5.0, 3.0, 6.0, 5.0, 2.0, 2.0, 9.0, 5.0, 3.0, 7.0, 3.0, 4.0, 2.0, 10.0, 1.0, 7.0, 2.0, 7.0, 1.0, 6.0, 8.0, 10.0, 1.0, 3.0, 3.0, 8.0, 5.0, 7.0, 10.0, 8.0, 3.0, 4.0, 4.0, 6.0, 2.0, 4.0, 2.0, 8.0, 9.0, 7.0, 2.0, 1.0, 7.0, 8.0, 9.0, 3.0, 9.0, 3.0, 5.0, 10.0, 1.0, 3.0, 7.0, 3.0, 5.0, 8.0, 1.0, 7.0, 10.0, 7.0, 8.0, 6.0, 10.0, 6.0, 3.0, 7.0, 6.0, 4.0, 7.0, 5.0, 10.0, 10.0, 10.0, 10.0, 6.0, 3.0, 10.0, 9.0, 6.0, 2.0, 5.0, 6.0, 3.0, 8.0, 5.0, 6.0, 2.0, 10.0, 3.0, 10.0, 1.0, 7.0, 2.0, 2.0, 5.0, 7.0, 5.0, 4.0, 10.0, 5.0, 2.0, 8.0, 9.0, 7.0, 9.0, 6.0, 10.0, 9.0, 6.0, 4.0, 6.0, 7.0, 6.0, 4.0, 10.0, 6.0, 8.0, 6.0, 4.0, 3.0, 3.0, 6.0, 10.0, 2.0, 10.0, 3.0, 8.0, 4.0, 6.0, 1.0, 7.0, 2.0, 6.0, 8.0, 2.0, 8.0, 7.0, 1.0, 7.0, 6.0, 6.0, 10.0, 7.0, 5.0, 2.0, 1.0, 5.0, 2.0, 8.0, 2.0, 3.0, 6.0, 5.0, 2.0, 10.0, 5.0, 4.0, 10.0, 6.0, 3.0, 4.0, 7.0, 8.0, 9.0, 2.0, 5.0, 2.0, 2.0, 9.0, 5.0, 8.0, 8.0, 5.0, 2.0, 2.0, 4.0, 6.0, 8.0, 2.0, 4.0, 7.0, 9.0, 5.0, 8.0, 4.0, 9.0, 1.0, 3.0, 3.0, 2.0, 7.0, 7.0, 7.0, 1.0, 1.0, 1.0, 5.0, 5.0, 7.0, 3.0, 2.0, 6.0, 9.0, 2.0, 2.0, 3.0, 7.0]
global b_x = 5
global d_y = [4.0, 3.0, 1.0, 6.0, 7.0, 5.0, 6.0, 8.0, 10.0, 2.0, 5.0, 8.0, 1.0, 7.0, 4.0, 2.0, 2.0, 2.0, 10.0, 5.0, 8.0, 2.0, 3.0, 8.0, 8.0, 2.0, 8.0, 2.0, 5.0, 6.0, 6.0, 4.0, 10.0, 3.0, 1.0, 6.0, 6.0, 7.0, 1.0, 5.0, 3.0, 4.0, 4.0, 8.0, 2.0, 8.0, 2.0, 3.0, 6.0, 3.0, 1.0, 9.0, 5.0, 1.0, 6.0, 8.0, 10.0, 5.0, 10.0, 6.0, 10.0, 2.0, 9.0, 4.0, 3.0, 5.0, 10.0, 7.0, 3.0, 2.0, 10.0, 7.0, 7.0, 1.0, 8.0, 1.0, 10.0, 10.0, 4.0, 10.0, 2.0, 6.0, 10.0, 10.0, 7.0, 4.0, 5.0, 9.0, 1.0, 10.0, 2.0, 9.0, 9.0, 7.0, 5.0, 3.0, 3.0, 5.0, 5.0, 6.0, 7.0, 9.0, 5.0, 4.0, 1.0, 7.0, 4.0, 9.0, 6.0, 9.0, 8.0, 3.0, 3.0, 10.0, 4.0, 5.0, 9.0, 7.0, 1.0, 1.0, 10.0, 8.0, 3.0, 2.0, 8.0, 9.0, 7.0, 8.0, 8.0, 10.0, 8.0, 5.0, 3.0, 5.0, 7.0, 2.0, 8.0, 10.0, 4.0, 9.0, 8.0, 8.0, 3.0, 9.0, 10.0, 7.0, 9.0, 8.0, 7.0, 2.0, 9.0, 1.0, 7.0, 9.0, 3.0, 5.0, 9.0, 3.0, 2.0, 3.0, 1.0, 8.0, 1.0, 3.0, 2.0, 8.0, 2.0, 7.0, 6.0, 3.0, 3.0, 7.0, 4.0, 4.0, 6.0, 6.0, 10.0, 2.0, 9.0, 5.0, 4.0, 6.0, 2.0, 1.0, 3.0, 3.0, 5.0, 3.0, 10.0, 4.0, 4.0, 10.0, 8.0, 9.0, 7.0, 4.0, 1.0, 2.0, 2.0, 3.0, 9.0, 3.0, 4.0, 8.0, 7.0, 1.0, 10.0, 5.0, 10.0, 5.0, 8.0, 3.0, 10.0, 10.0, 9.0, 8.0, 2.0, 4.0, 4.0, 6.0, 1.0, 4.0, 9.0, 5.0, 5.0, 7.0, 10.0, 10.0, 3.0, 1.0, 8.0, 5.0, 8.0, 8.0, 5.0, 3.0, 2.0, 3.0, 3.0, 1.0, 9.0, 4.0, 6.0, 8.0, 1.0, 10.0, 6.0, 5.0, 6.0, 1.0, 4.0, 6.0, 5.0, 2.0, 1.0, 9.0, 1.0, 7.0, 9.0, 5.0, 5.0, 1.0, 5.0, 2.0, 6.0, 5.0, 3.0, 7.0, 4.0, 5.0, 1.0, 5.0, 2.0, 7.0, 10.0, 8.0, 10.0, 2.0, 7.0, 1.0, 7.0, 7.0, 9.0, 9.0, 8.0, 4.0, 8.0, 10.0, 2.0, 5.0, 9.0, 5.0, 3.0, 5.0, 10.0, 9.0, 1.0, 7.0, 1.0, 1.0, 6.0, 10.0, 7.0, 3.0, 6.0, 10.0, 4.0, 8.0, 2.0, 8.0, 9.0, 10.0, 9.0, 8.0, 9.0, 1.0, 2.0, 4.0, 10.0, 9.0, 1.0, 7.0, 4.0, 8.0, 6.0, 1.0, 9.0, 5.0, 5.0, 8.0, 2.0, 3.0, 4.0, 10.0, 10.0, 2.0, 6.0, 6.0, 3.0, 6.0, 7.0, 8.0]
global b_y = 10
global p = [0.113, 0.573, 0.499, 0.962, 0.672, 0.781, 0.011, 0.82, 0.931, 0.352, 0.595, 0.088, 0.54, 0.608, 0.159, 0.968, 0.011, 0.822, 0.795, 0.332, 0.036, 0.225, 0.207, 0.453, 0.615, 0.952, 0.106, 0.288, 0.767, 0.289, 0.85, 0.306, 0.154, 0.012, 0.734, 0.113, 0.62, 0.85, 0.329, 0.19, 0.511, 0.672, 0.714, 0.307, 0.058, 0.633, 0.752, 0.557, 0.135, 0.023, 0.654, 0.102, 0.08, 0.741, 0.78, 0.498, 0.441, 0.777, 0.358, 0.933, 0.976, 0.489, 0.371, 0.022, 0.364, 0.375, 0.259, 0.611, 0.295, 0.855, 0.548, 0.017, 0.195, 0.096, 0.263, 0.426, 0.152, 0.134, 0.185, 0.711, 0.707, 0.452, 0.858, 0.652, 0.685, 0.205, 0.997, 0.903, 0.545, 0.341, 0.058, 0.017, 0.138, 0.052, 0.58, 0.304, 0.184, 0.814, 0.571, 0.521, 0.001, 0.122, 0.672, 0.596, 0.825, 0.138, 0.547, 0.409, 0.478, 0.344, 0.58, 0.013, 0.112, 0.201, 0.283, 0.887, 0.883, 0.415, 0.834, 0.038, 0.218, 0.367, 0.584, 0.281, 0.636, 0.319, 0.964, 0.593, 0.865, 0.437, 0.632, 0.368, 0.768, 0.591, 0.766, 0.24, 0.601, 0.374, 0.298, 0.767, 0.807, 0.849, 0.968, 0.256, 0.114, 0.825, 0.805, 0.306, 0.638, 0.08, 0.854, 0.418, 0.684, 0.05, 0.742, 0.503, 0.265, 0.719, 0.133, 0.458, 0.114, 0.369, 0.583, 0.84, 0.298, 0.144, 0.562, 0.371, 0.287, 0.423, 0.856, 0.16, 0.72, 0.925, 0.028, 0.918, 0.61, 0.456, 0.152, 0.028, 0.449, 0.847, 0.664, 0.607, 0.736, 0.149, 0.27, 0.961, 0.157, 0.321, 0.304, 0.326, 0.632, 0.45, 0.539, 0.122, 0.977, 0.379, 0.285, 0.027, 0.141, 0.469, 0.187, 0.39, 0.598, 0.574, 0.147, 0.265, 0.457, 0.893, 0.117, 0.121, 0.224, 0.448, 0.382, 0.456, 0.778, 0.406, 0.598, 0.779, 0.968, 0.118, 0.712, 0.762, 0.92, 0.256, 0.299, 0.119, 0.173, 0.6, 0.889, 0.766, 0.295, 0.163, 0.371, 0.564, 0.702, 0.575, 0.952, 0.909, 0.963, 0.285, 0.149, 0.247, 0.969, 0.921, 0.861, 0.556, 0.565, 0.134, 0.9, 0.557, 0.609, 0.646, 0.222, 0.229, 0.173, 0.25, 0.007, 0.056, 0.646, 0.843, 0.971, 0.147, 0.539, 0.582, 0.922, 0.705, 0.214, 0.127, 0.47, 0.897, 0.235, 0.273, 0.84, 0.252, 0.718, 0.607, 0.348, 0.596, 0.533, 0.873, 0.273, 0.602, 0.084, 0.108, 0.505, 0.528, 0.562, 0.097, 0.701, 0.587, 0.343, 0.912, 0.649, 0.921, 0.729, 0.614, 0.418, 0.488, 0.719, 0.277, 0.662, 0.992, 0.613, 0.495, 0.05, 0.986, 0.196, 0.001, 0.475, 0.35, 0.9, 0.142, 0.959, 0.544, 0.002, 0.439, 0.257, 0.847, 0.257, 0.234, 0.618, 0.589, 0.529, 0.524, 0.731, 0.328, 0.863, 0.272, 0.667, 0.853, 0.848, 0.709, 0.913, 0.766, 0.3, 0.556, 0.822, 0.272, 0.671, 0.211]
global q = [0.716, 0.717, 0.556, 0.965, 0.887, 0.891, 0.561, 0.991, 0.963, 0.635, 0.995, 0.333, 0.659, 0.817, 0.629, 0.985, 0.574, 0.858, 0.849, 0.621, 0.495, 0.923, 0.491, 0.499, 0.75, 0.961, 0.391, 0.574, 0.851, 0.828, 0.938, 0.804, 0.673, 0.638, 0.775, 0.926, 0.661, 0.939, 0.486, 0.225, 0.938, 0.836, 0.825, 0.923, 0.333, 0.657, 0.939, 0.589, 0.884, 0.94, 0.737, 0.696, 0.94, 0.804, 0.871, 0.819, 0.947, 0.865, 0.463, 0.981, 0.99, 0.513, 0.898, 0.552, 0.463, 0.769, 0.477, 0.719, 0.469, 0.937, 0.566, 0.464, 0.427, 0.797, 0.572, 0.828, 0.932, 0.624, 0.325, 0.908, 0.903, 0.698, 0.901, 0.797, 0.694, 0.368, 0.997, 0.938, 0.574, 0.711, 0.856, 0.759, 0.161, 0.601, 0.701, 0.929, 0.563, 0.893, 0.691, 0.82, 0.734, 0.719, 0.715, 0.787, 0.932, 0.272, 0.636, 0.688, 0.967, 0.968, 0.625, 0.625, 0.37, 0.829, 0.375, 0.908, 0.95, 0.958, 0.859, 0.649, 0.899, 0.495, 0.584, 0.832, 0.734, 0.935, 0.965, 0.833, 0.987, 0.907, 0.957, 0.6, 0.979, 0.688, 0.824, 0.534, 0.899, 0.539, 0.459, 0.88, 0.964, 0.947, 0.976, 0.843, 0.793, 0.856, 0.926, 0.479, 0.892, 0.466, 0.941, 0.879, 0.896, 0.444, 0.991, 0.505, 0.531, 0.735, 0.909, 0.558, 0.606, 0.713, 0.928, 0.928, 0.602, 0.218, 0.812, 0.593, 0.355, 0.451, 0.943, 0.586, 0.906, 0.968, 0.122, 0.931, 0.962, 0.565, 0.487, 0.975, 0.912, 0.948, 0.665, 0.988, 0.983, 0.513, 0.88, 0.979, 0.601, 0.527, 0.724, 0.749, 0.756, 0.823, 0.661, 0.968, 0.981, 0.89, 0.416, 0.485, 0.868, 0.886, 0.572, 0.673, 0.615, 0.683, 0.685, 0.716, 0.953, 0.918, 0.306, 0.663, 0.782, 0.724, 0.734, 0.559, 0.788, 0.937, 0.68, 0.844, 0.997, 0.691, 0.785, 0.805, 0.973, 0.479, 0.867, 0.267, 0.179, 0.856, 0.96, 0.932, 0.352, 0.773, 0.894, 0.886, 0.945, 0.685, 0.978, 0.964, 0.978, 0.667, 0.864, 0.454, 0.97, 0.926, 0.88, 0.817, 0.774, 0.912, 0.953, 0.996, 0.745, 0.697, 0.986, 0.687, 0.456, 0.357, 0.139, 0.545, 0.806, 0.934, 0.985, 0.858, 0.554, 0.613, 0.922, 0.981, 0.485, 0.223, 0.832, 0.939, 0.545, 0.765, 0.929, 0.329, 0.766, 0.884, 0.547, 0.702, 0.984, 0.894, 0.434, 0.957, 0.434, 0.988, 0.732, 0.573, 0.746, 0.281, 0.813, 0.959, 0.446, 0.96, 0.784, 0.989, 0.893, 0.662, 0.54, 0.98, 0.968, 0.464, 0.974, 0.996, 0.739, 0.702, 0.849, 0.999, 0.802, 0.852, 0.699, 0.507, 0.995, 0.944, 0.963, 0.71, 0.767, 0.504, 0.455, 0.938, 0.552, 0.499, 0.717, 0.812, 0.759, 0.895, 0.964, 0.585, 0.952, 0.316, 0.816, 0.977, 0.896, 0.891, 0.99, 0.809, 0.75, 0.799, 0.935, 0.876, 0.72, 0.57]
global origin = 1
global destination = 60