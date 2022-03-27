global arcs = [1 4; 1 6; 1 9; 1 39; 1 48; 2 11; 2 13; 2 15; 2 27; 2 30; 2 43; 2 49; 3 8; 3 38; 3 40; 3 50; 4 14; 4 46; 5 3; 5 15; 5 28; 5 32; 5 34; 5 44; 5 45; 6 2; 6 25; 6 28; 6 39; 7 2; 7 13; 7 20; 7 28; 7 33; 7 46; 8 36; 8 42; 9 13; 9 21; 9 27; 10 2; 10 8; 10 17; 10 23; 10 35; 10 38; 10 39; 10 46; 10 49; 11 15; 11 36; 11 42; 12 13; 12 24; 12 47; 13 5; 13 8; 13 12; 13 27; 13 31; 13 33; 13 39; 13 47; 14 11; 14 12; 14 17; 14 35; 14 49; 15 17; 15 20; 15 27; 15 32; 15 34; 15 43; 16 5; 16 11; 16 27; 17 10; 17 19; 17 24; 17 45; 17 47; 18 19; 18 45; 19 25; 19 30; 19 44; 20 11; 20 31; 20 33; 20 34; 20 37; 20 38; 20 44; 20 49; 21 7; 21 12; 21 14; 21 22; 21 33; 21 48; 22 4; 22 7; 22 9; 22 11; 22 13; 22 27; 22 29; 22 30; 22 31; 22 45; 23 13; 23 49; 24 23; 24 25; 24 29; 25 5; 25 7; 25 14; 25 15; 25 50; 26 8; 26 31; 26 40; 27 11; 27 17; 27 28; 28 23; 28 36; 29 3; 29 11; 29 22; 29 24; 29 35; 30 11; 30 24; 30 35; 30 42; 30 44; 30 50; 31 17; 31 24; 31 30; 31 44; 32 19; 32 38; 32 42; 33 6; 33 9; 33 18; 33 23; 33 24; 33 25; 33 29; 34 15; 34 29; 34 42; 34 47; 34 49; 34 50; 35 28; 36 11; 36 26; 36 48; 37 12; 37 18; 37 25; 38 16; 38 20; 38 32; 38 43; 38 46; 39 3; 39 4; 39 7; 39 8; 39 12; 39 42; 40 6; 40 13; 40 25; 40 28; 40 38; 41 5; 41 19; 42 5; 42 16; 42 36; 42 39; 43 4; 43 17; 43 22; 43 37; 43 47; 43 48; 44 6; 44 7; 44 24; 44 34; 44 37; 44 49; 45 24; 46 12; 46 35; 46 39; 47 25; 47 26; 47 34; 47 38; 47 45; 47 46; 48 3; 48 6; 48 26; 48 30; 48 36; 48 46; 49 5; 49 29; 49 39; 49 41; 49 48]
global d_x = [10.0, 10.0, 6.0, 2.0, 4.0, 10.0, 7.0, 4.0, 2.0, 2.0, 6.0, 2.0, 2.0, 4.0, 10.0, 4.0, 3.0, 2.0, 6.0, 7.0, 5.0, 2.0, 10.0, 7.0, 7.0, 6.0, 4.0, 6.0, 4.0, 9.0, 6.0, 1.0, 1.0, 1.0, 2.0, 9.0, 9.0, 8.0, 4.0, 7.0, 4.0, 6.0, 4.0, 7.0, 8.0, 7.0, 8.0, 8.0, 7.0, 8.0, 2.0, 3.0, 9.0, 5.0, 1.0, 6.0, 5.0, 6.0, 8.0, 1.0, 3.0, 7.0, 1.0, 10.0, 2.0, 1.0, 7.0, 1.0, 4.0, 2.0, 3.0, 2.0, 3.0, 6.0, 8.0, 8.0, 4.0, 7.0, 6.0, 6.0, 4.0, 6.0, 4.0, 3.0, 2.0, 10.0, 7.0, 9.0, 4.0, 1.0, 1.0, 5.0, 6.0, 10.0, 3.0, 6.0, 4.0, 8.0, 10.0, 8.0, 1.0, 7.0, 4.0, 6.0, 9.0, 7.0, 4.0, 9.0, 5.0, 1.0, 9.0, 3.0, 8.0, 5.0, 2.0, 5.0, 6.0, 2.0, 9.0, 8.0, 3.0, 2.0, 8.0, 3.0, 3.0, 3.0, 7.0, 9.0, 9.0, 4.0, 5.0, 9.0, 5.0, 10.0, 1.0, 3.0, 4.0, 8.0, 5.0, 7.0, 7.0, 2.0, 2.0, 8.0, 6.0, 10.0, 1.0, 6.0, 7.0, 8.0, 6.0, 3.0, 10.0, 3.0, 8.0, 1.0, 3.0, 1.0, 5.0, 10.0, 5.0, 9.0, 7.0, 10.0, 9.0, 10.0, 3.0, 10.0, 2.0, 2.0, 9.0, 8.0, 7.0, 9.0, 7.0, 1.0, 9.0, 6.0, 2.0, 8.0, 3.0, 1.0, 5.0, 1.0, 9.0, 3.0, 3.0, 1.0, 7.0, 5.0, 2.0, 4.0, 1.0, 2.0, 3.0, 9.0, 2.0, 4.0, 1.0, 6.0, 6.0, 2.0, 8.0, 2.0, 6.0, 8.0, 5.0, 2.0, 5.0, 4.0, 1.0, 10.0, 5.0, 4.0, 4.0, 9.0, 5.0, 3.0, 1.0, 1.0, 1.0, 2.0]
global b_x = 5
global d_y = [3.0, 8.0, 9.0, 2.0, 7.0, 6.0, 9.0, 3.0, 2.0, 9.0, 5.0, 5.0, 8.0, 9.0, 4.0, 10.0, 10.0, 6.0, 10.0, 10.0, 1.0, 2.0, 5.0, 4.0, 2.0, 6.0, 4.0, 1.0, 4.0, 2.0, 10.0, 4.0, 5.0, 6.0, 5.0, 1.0, 8.0, 9.0, 10.0, 1.0, 9.0, 3.0, 1.0, 5.0, 8.0, 10.0, 8.0, 1.0, 4.0, 9.0, 7.0, 4.0, 9.0, 5.0, 3.0, 1.0, 8.0, 10.0, 1.0, 10.0, 7.0, 6.0, 3.0, 8.0, 6.0, 5.0, 8.0, 6.0, 7.0, 8.0, 4.0, 8.0, 3.0, 1.0, 7.0, 7.0, 4.0, 5.0, 10.0, 2.0, 6.0, 9.0, 3.0, 5.0, 7.0, 6.0, 10.0, 7.0, 7.0, 4.0, 10.0, 8.0, 5.0, 7.0, 2.0, 8.0, 3.0, 9.0, 1.0, 10.0, 2.0, 7.0, 9.0, 8.0, 7.0, 2.0, 2.0, 2.0, 1.0, 10.0, 4.0, 10.0, 7.0, 7.0, 2.0, 10.0, 10.0, 8.0, 4.0, 8.0, 4.0, 9.0, 7.0, 2.0, 1.0, 6.0, 3.0, 2.0, 6.0, 9.0, 4.0, 8.0, 3.0, 3.0, 4.0, 1.0, 3.0, 8.0, 5.0, 5.0, 5.0, 1.0, 7.0, 3.0, 7.0, 5.0, 3.0, 9.0, 10.0, 1.0, 9.0, 7.0, 10.0, 10.0, 7.0, 3.0, 9.0, 8.0, 10.0, 9.0, 9.0, 2.0, 8.0, 5.0, 2.0, 2.0, 8.0, 3.0, 6.0, 7.0, 6.0, 2.0, 1.0, 5.0, 5.0, 3.0, 6.0, 9.0, 8.0, 6.0, 6.0, 10.0, 5.0, 7.0, 6.0, 10.0, 8.0, 6.0, 7.0, 4.0, 5.0, 6.0, 10.0, 8.0, 1.0, 2.0, 5.0, 5.0, 2.0, 5.0, 3.0, 9.0, 2.0, 6.0, 1.0, 4.0, 9.0, 10.0, 9.0, 7.0, 8.0, 4.0, 10.0, 5.0, 1.0, 4.0, 1.0, 10.0, 7.0, 8.0, 4.0, 9.0]
global b_y = 10
global p = [0.69, 0.525, 0.14, 0.135, 0.707, 0.223, 0.426, 0.948, 0.793, 0.846, 0.73, 0.656, 0.572, 0.923, 0.705, 0.921, 0.775, 0.612, 0.036, 0.973, 0.914, 0.468, 0.301, 0.268, 0.071, 0.102, 0.803, 0.391, 0.033, 0.044, 0.066, 0.14, 0.445, 0.91, 0.063, 0.746, 0.835, 0.001, 0.326, 0.066, 0.184, 0.659, 0.537, 0.757, 0.937, 0.85, 0.247, 0.148, 0.81, 0.564, 0.483, 0.122, 0.234, 0.578, 0.299, 0.497, 0.575, 0.847, 0.901, 0.118, 0.788, 0.22, 0.763, 0.05, 0.644, 0.543, 0.731, 0.689, 0.797, 0.751, 0.34, 0.499, 0.929, 0.798, 0.275, 0.272, 0.954, 0.671, 0.944, 0.061, 0.709, 0.691, 0.877, 0.684, 0.527, 0.364, 0.705, 0.003, 0.339, 0.789, 0.265, 0.04, 0.708, 0.813, 0.319, 0.264, 0.319, 0.014, 0.036, 0.867, 0.942, 0.749, 0.246, 0.993, 0.904, 0.766, 0.563, 0.178, 0.642, 0.143, 0.615, 0.156, 0.134, 0.426, 0.333, 0.822, 0.761, 0.05, 0.182, 0.849, 0.85, 0.479, 0.938, 0.881, 0.558, 0.581, 0.414, 0.032, 0.715, 0.187, 0.53, 0.411, 0.119, 0.723, 0.027, 0.89, 0.763, 0.268, 0.964, 0.795, 0.148, 0.462, 0.86, 0.255, 0.399, 0.932, 0.512, 0.428, 0.445, 0.666, 0.06, 0.541, 0.289, 0.536, 0.875, 0.586, 0.659, 0.569, 0.526, 0.033, 0.941, 0.403, 0.648, 0.789, 0.137, 0.356, 0.307, 0.684, 0.472, 0.8, 0.514, 0.309, 0.417, 0.313, 0.42, 0.657, 0.94, 0.428, 0.276, 0.417, 0.587, 0.82, 0.481, 0.829, 0.488, 0.851, 0.262, 0.986, 0.401, 0.667, 0.573, 0.276, 0.734, 0.358, 0.731, 0.483, 0.01, 0.521, 0.193, 0.594, 0.673, 0.7, 0.207, 0.792, 0.34, 0.487, 0.889, 0.107, 0.286, 0.13, 0.224, 0.428, 0.901, 0.296, 0.726, 0.74, 0.116, 0.73, 0.141, 0.877, 0.798, 0.171]
global q = [0.898, 0.874, 0.597, 0.407, 0.918, 0.486, 0.602, 0.993, 0.902, 0.937, 0.858, 0.851, 0.915, 0.991, 0.906, 0.977, 0.822, 0.874, 0.245, 0.979, 0.988, 0.842, 0.829, 0.387, 0.836, 0.784, 0.96, 0.524, 0.209, 0.096, 0.532, 0.381, 0.455, 0.967, 0.721, 0.957, 0.923, 0.221, 0.535, 0.37, 0.509, 0.861, 0.889, 0.908, 0.956, 0.895, 0.813, 0.419, 0.824, 0.747, 0.636, 0.646, 0.489, 0.812, 0.908, 0.694, 0.888, 0.87, 0.914, 0.171, 0.819, 0.257, 0.8, 0.822, 0.75, 0.963, 0.896, 0.726, 0.919, 0.951, 0.474, 0.734, 0.934, 0.95, 0.355, 0.474, 0.973, 0.767, 0.992, 0.708, 0.757, 0.873, 0.903, 0.761, 0.695, 0.521, 0.844, 0.444, 0.816, 0.965, 0.636, 0.445, 0.762, 0.831, 0.732, 0.54, 0.712, 0.06, 0.3, 0.87, 0.975, 0.964, 0.75, 0.999, 0.92, 0.989, 0.571, 0.455, 0.97, 0.397, 0.684, 0.538, 0.876, 0.864, 0.689, 0.922, 0.918, 0.7, 0.969, 0.861, 0.993, 0.919, 0.941, 0.906, 0.724, 0.634, 0.572, 0.439, 0.897, 0.941, 0.578, 0.953, 0.148, 0.952, 0.551, 0.951, 0.946, 0.848, 0.966, 0.885, 0.66, 0.914, 0.881, 0.923, 0.811, 0.979, 0.966, 0.457, 0.766, 0.676, 0.849, 0.807, 0.652, 0.835, 0.897, 0.861, 0.731, 0.906, 0.827, 0.685, 0.944, 0.456, 0.976, 0.894, 0.546, 0.433, 0.756, 0.706, 0.933, 0.857, 0.876, 0.819, 0.66, 0.747, 0.677, 0.78, 0.99, 0.901, 0.638, 0.969, 0.907, 0.855, 0.522, 0.891, 0.683, 0.967, 0.547, 0.997, 0.788, 0.87, 0.942, 0.707, 0.951, 0.389, 0.952, 0.83, 0.895, 0.954, 0.293, 0.694, 0.992, 0.796, 0.883, 0.93, 0.38, 0.566, 0.95, 0.96, 0.766, 0.339, 0.319, 0.661, 0.938, 0.871, 0.733, 0.786, 0.65, 0.824, 0.435, 0.899, 0.894, 0.443]
global origin = 1
global destination = 50