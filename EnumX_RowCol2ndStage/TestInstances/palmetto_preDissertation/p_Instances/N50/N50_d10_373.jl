global arcs = [1 3; 1 6; 1 11; 1 23; 2 8; 2 13; 2 22; 2 35; 2 41; 3 6; 3 22; 3 25; 3 35; 4 25; 4 28; 5 10; 5 12; 5 21; 5 40; 5 46; 6 2; 6 5; 6 11; 6 27; 6 32; 6 40; 7 2; 7 5; 7 8; 7 33; 7 34; 7 36; 7 39; 8 6; 8 11; 8 18; 8 20; 8 21; 9 6; 9 23; 9 27; 9 42; 9 50; 10 12; 10 23; 10 32; 10 36; 11 3; 11 28; 11 30; 11 36; 11 40; 12 8; 12 20; 12 26; 12 33; 12 37; 12 43; 12 44; 13 9; 13 15; 13 26; 13 27; 13 29; 13 30; 14 3; 14 19; 14 33; 14 47; 15 12; 15 13; 15 24; 15 25; 15 33; 15 35; 15 36; 15 39; 15 41; 15 49; 16 2; 16 6; 16 20; 16 33; 16 41; 16 42; 16 43; 17 2; 17 3; 17 6; 17 26; 17 29; 17 43; 18 12; 18 13; 18 37; 18 42; 18 46; 19 7; 19 27; 19 42; 19 50; 20 6; 20 10; 20 15; 20 16; 20 17; 20 27; 20 39; 20 40; 20 42; 20 46; 20 48; 21 14; 21 17; 21 30; 21 39; 21 41; 21 47; 22 5; 22 14; 22 26; 22 38; 22 48; 23 8; 23 27; 23 34; 23 44; 24 16; 24 31; 24 38; 24 41; 24 49; 25 6; 25 9; 25 26; 25 28; 25 31; 25 35; 25 39; 25 43; 25 44; 25 50; 26 14; 26 21; 26 23; 26 37; 26 39; 27 14; 27 30; 27 31; 28 6; 28 7; 28 8; 28 19; 28 23; 28 38; 28 42; 29 15; 29 17; 29 25; 29 39; 29 45; 29 46; 30 6; 30 10; 30 14; 30 15; 30 26; 30 48; 31 12; 31 13; 31 39; 31 41; 32 8; 32 13; 32 15; 32 19; 32 31; 32 36; 33 8; 33 14; 33 35; 33 45; 34 18; 34 19; 34 37; 35 25; 35 26; 35 42; 35 46; 36 23; 36 28; 36 31; 36 35; 36 39; 36 44; 36 47; 37 6; 37 27; 37 44; 37 49; 38 2; 38 4; 38 10; 38 48; 39 3; 39 8; 39 10; 39 12; 39 24; 39 27; 39 37; 39 41; 39 44; 40 3; 40 8; 40 19; 40 35; 41 10; 41 22; 41 33; 41 46; 42 4; 42 20; 42 22; 42 25; 42 36; 42 44; 42 45; 42 48; 43 3; 43 10; 43 11; 43 15; 43 18; 43 23; 43 30; 43 46; 44 12; 44 13; 44 15; 44 32; 44 33; 44 34; 44 38; 45 16; 45 20; 45 37; 46 13; 46 17; 46 22; 46 34; 46 38; 46 47; 47 12; 47 19; 47 32; 48 6; 48 18; 48 19; 48 24; 48 25; 48 28; 48 31; 48 40; 49 3; 49 28; 49 33]
global d_x = [4.0, 4.0, 10.0, 7.0, 7.0, 10.0, 6.0, 1.0, 9.0, 10.0, 10.0, 3.0, 8.0, 1.0, 2.0, 6.0, 10.0, 4.0, 2.0, 2.0, 5.0, 10.0, 1.0, 9.0, 8.0, 8.0, 2.0, 1.0, 8.0, 3.0, 9.0, 10.0, 1.0, 6.0, 9.0, 3.0, 3.0, 4.0, 10.0, 4.0, 8.0, 2.0, 9.0, 5.0, 5.0, 8.0, 2.0, 1.0, 7.0, 5.0, 6.0, 10.0, 2.0, 6.0, 7.0, 6.0, 1.0, 1.0, 6.0, 5.0, 3.0, 6.0, 7.0, 5.0, 10.0, 7.0, 2.0, 10.0, 5.0, 6.0, 6.0, 9.0, 10.0, 8.0, 4.0, 4.0, 7.0, 3.0, 6.0, 7.0, 4.0, 1.0, 2.0, 8.0, 9.0, 6.0, 8.0, 3.0, 5.0, 1.0, 2.0, 10.0, 4.0, 6.0, 6.0, 6.0, 3.0, 5.0, 10.0, 6.0, 3.0, 2.0, 4.0, 9.0, 9.0, 10.0, 2.0, 9.0, 6.0, 4.0, 8.0, 8.0, 5.0, 4.0, 9.0, 9.0, 8.0, 2.0, 8.0, 1.0, 3.0, 10.0, 5.0, 5.0, 10.0, 2.0, 10.0, 7.0, 6.0, 10.0, 10.0, 2.0, 3.0, 9.0, 6.0, 10.0, 6.0, 10.0, 7.0, 1.0, 5.0, 7.0, 9.0, 1.0, 6.0, 7.0, 4.0, 1.0, 4.0, 4.0, 4.0, 10.0, 8.0, 3.0, 7.0, 7.0, 4.0, 5.0, 7.0, 1.0, 6.0, 4.0, 8.0, 7.0, 5.0, 10.0, 3.0, 5.0, 8.0, 8.0, 8.0, 5.0, 1.0, 5.0, 9.0, 7.0, 1.0, 2.0, 6.0, 5.0, 4.0, 2.0, 3.0, 10.0, 4.0, 5.0, 4.0, 7.0, 1.0, 4.0, 6.0, 3.0, 9.0, 6.0, 4.0, 2.0, 6.0, 8.0, 6.0, 4.0, 9.0, 5.0, 9.0, 6.0, 8.0, 10.0, 5.0, 8.0, 1.0, 6.0, 4.0, 3.0, 3.0, 6.0, 3.0, 5.0, 5.0, 9.0, 3.0, 6.0, 5.0, 3.0, 10.0, 5.0, 5.0, 3.0, 2.0, 9.0, 6.0, 3.0, 2.0, 6.0, 4.0, 2.0, 2.0, 8.0, 9.0, 1.0, 8.0, 4.0, 8.0, 1.0, 4.0, 6.0, 3.0, 3.0, 1.0, 8.0, 2.0, 10.0, 9.0, 10.0, 10.0, 10.0, 6.0, 9.0, 4.0, 5.0, 10.0, 6.0, 6.0, 1.0, 8.0, 3.0, 3.0, 4.0, 5.0, 7.0]
global b_x = 5
global d_y = [6.0, 8.0, 6.0, 5.0, 6.0, 8.0, 2.0, 9.0, 5.0, 7.0, 6.0, 9.0, 4.0, 10.0, 8.0, 1.0, 3.0, 7.0, 7.0, 6.0, 7.0, 2.0, 1.0, 1.0, 9.0, 6.0, 1.0, 2.0, 9.0, 4.0, 9.0, 5.0, 4.0, 8.0, 9.0, 9.0, 7.0, 6.0, 4.0, 10.0, 3.0, 1.0, 3.0, 2.0, 9.0, 7.0, 2.0, 10.0, 10.0, 4.0, 4.0, 9.0, 5.0, 5.0, 7.0, 7.0, 8.0, 3.0, 9.0, 5.0, 6.0, 6.0, 8.0, 8.0, 2.0, 8.0, 1.0, 3.0, 8.0, 6.0, 8.0, 2.0, 1.0, 10.0, 6.0, 8.0, 3.0, 8.0, 5.0, 6.0, 2.0, 6.0, 8.0, 10.0, 1.0, 10.0, 8.0, 10.0, 7.0, 5.0, 10.0, 4.0, 4.0, 9.0, 6.0, 8.0, 4.0, 5.0, 1.0, 4.0, 2.0, 8.0, 1.0, 3.0, 10.0, 6.0, 2.0, 1.0, 5.0, 7.0, 7.0, 7.0, 1.0, 10.0, 1.0, 7.0, 6.0, 9.0, 7.0, 3.0, 2.0, 4.0, 3.0, 5.0, 7.0, 5.0, 2.0, 8.0, 8.0, 6.0, 6.0, 3.0, 6.0, 4.0, 5.0, 5.0, 8.0, 4.0, 9.0, 10.0, 7.0, 7.0, 10.0, 2.0, 4.0, 1.0, 5.0, 9.0, 1.0, 10.0, 8.0, 2.0, 8.0, 6.0, 5.0, 3.0, 1.0, 4.0, 6.0, 1.0, 5.0, 9.0, 3.0, 8.0, 8.0, 10.0, 1.0, 5.0, 10.0, 1.0, 5.0, 1.0, 5.0, 1.0, 4.0, 4.0, 8.0, 7.0, 8.0, 2.0, 10.0, 10.0, 7.0, 6.0, 1.0, 4.0, 4.0, 10.0, 6.0, 5.0, 2.0, 2.0, 7.0, 4.0, 5.0, 6.0, 2.0, 6.0, 6.0, 3.0, 7.0, 10.0, 6.0, 3.0, 1.0, 9.0, 3.0, 9.0, 6.0, 2.0, 8.0, 1.0, 10.0, 7.0, 1.0, 2.0, 10.0, 1.0, 10.0, 1.0, 4.0, 5.0, 4.0, 8.0, 9.0, 6.0, 7.0, 9.0, 9.0, 2.0, 9.0, 9.0, 3.0, 10.0, 4.0, 5.0, 6.0, 5.0, 3.0, 2.0, 4.0, 5.0, 10.0, 7.0, 5.0, 5.0, 1.0, 7.0, 5.0, 9.0, 6.0, 6.0, 5.0, 7.0, 7.0, 8.0, 8.0, 3.0, 3.0, 1.0, 9.0, 8.0, 3.0, 10.0, 2.0, 7.0, 1.0, 4.0]
global b_y = 10
global p = [0.524, 0.022, 0.326, 0.169, 0.578, 0.036, 0.527, 0.391, 0.973, 0.754, 0.976, 0.317, 0.806, 0.348, 0.058, 0.229, 0.044, 0.74, 0.83, 0.975, 0.595, 0.974, 0.197, 0.057, 0.402, 0.508, 0.096, 0.503, 0.531, 0.398, 0.218, 0.978, 0.776, 0.541, 0.779, 0.147, 0.725, 0.876, 0.126, 0.243, 0.82, 0.559, 0.717, 0.047, 0.521, 0.683, 0.447, 0.843, 0.414, 0.74, 0.292, 0.955, 0.97, 0.014, 0.543, 0.149, 0.735, 0.436, 0.324, 0.405, 0.105, 0.362, 0.523, 0.481, 0.087, 0.657, 0.058, 0.672, 0.433, 0.666, 0.029, 0.422, 0.826, 0.381, 0.911, 0.749, 0.791, 0.232, 0.957, 0.965, 0.975, 0.09, 0.778, 0.583, 0.138, 0.394, 0.19, 0.618, 0.773, 0.263, 0.387, 0.571, 0.327, 0.834, 0.903, 0.917, 0.519, 0.817, 0.589, 0.711, 0.803, 0.186, 0.567, 0.639, 0.544, 0.24, 0.16, 0.99, 0.044, 0.626, 0.763, 0.669, 0.209, 0.777, 0.327, 0.889, 0.996, 0.586, 0.405, 0.389, 0.136, 0.399, 0.179, 0.91, 0.419, 0.886, 0.61, 0.303, 0.765, 0.943, 0.926, 0.337, 0.709, 0.688, 0.934, 0.07, 0.438, 0.117, 0.576, 0.254, 0.956, 0.297, 0.103, 0.244, 0.626, 0.8, 0.125, 0.559, 0.969, 0.615, 0.429, 0.867, 0.642, 0.677, 0.165, 0.343, 0.457, 0.056, 0.393, 0.505, 0.263, 0.258, 0.245, 0.17, 0.601, 0.057, 0.828, 0.391, 0.85, 0.101, 0.062, 0.064, 0.881, 0.011, 0.037, 0.042, 0.273, 0.222, 0.362, 0.158, 0.927, 0.729, 0.566, 0.233, 0.835, 0.13, 0.609, 0.334, 0.72, 0.537, 0.612, 0.864, 0.068, 0.197, 0.182, 0.934, 0.449, 0.395, 0.833, 0.004, 0.004, 0.201, 0.671, 0.366, 0.573, 0.945, 0.409, 0.558, 0.249, 0.705, 0.296, 0.668, 0.62, 0.722, 0.896, 0.086, 0.171, 0.655, 0.641, 0.749, 0.112, 0.507, 0.92, 0.022, 0.57, 0.888, 0.589, 0.495, 0.231, 0.817, 0.423, 0.695, 0.486, 0.278, 0.022, 0.417, 0.927, 0.346, 0.611, 0.105, 0.858, 0.424, 0.307, 0.861, 0.664, 0.581, 0.48, 0.033, 0.47, 0.428, 0.214, 0.412, 0.592, 0.455, 0.956, 0.276, 0.578, 0.357, 0.978, 0.495, 0.935, 0.572, 0.479, 0.086, 0.579, 0.448, 0.159, 0.114]
global q = [0.798, 0.187, 0.956, 0.286, 0.943, 0.849, 0.772, 0.766, 0.994, 0.789, 0.983, 0.494, 0.898, 0.404, 0.658, 0.464, 0.499, 0.816, 0.873, 0.982, 0.599, 0.999, 0.388, 0.901, 0.781, 0.618, 0.972, 0.568, 0.625, 0.461, 0.333, 0.99, 0.967, 0.596, 0.928, 0.973, 0.774, 0.949, 0.489, 0.337, 0.866, 0.908, 0.752, 0.662, 0.998, 0.688, 0.917, 0.867, 0.479, 0.966, 0.873, 0.966, 0.995, 0.514, 0.683, 0.159, 0.85, 0.477, 0.962, 0.492, 0.982, 0.983, 0.591, 0.622, 0.109, 0.857, 0.921, 0.741, 0.933, 0.98, 0.217, 0.957, 0.888, 0.646, 0.941, 0.976, 0.874, 0.557, 0.975, 0.989, 0.987, 0.76, 0.821, 0.842, 0.478, 0.926, 0.289, 0.97, 0.931, 0.319, 0.918, 0.674, 0.946, 0.912, 0.961, 0.949, 0.992, 0.963, 0.675, 0.95, 0.913, 0.222, 0.967, 0.865, 0.639, 0.292, 0.548, 0.998, 0.621, 0.635, 0.93, 0.968, 0.677, 0.924, 0.938, 0.958, 0.998, 0.714, 0.688, 0.933, 0.96, 0.848, 0.574, 0.993, 0.876, 0.9, 0.61, 0.536, 0.921, 0.951, 0.996, 0.747, 0.785, 0.912, 0.984, 0.462, 0.593, 0.866, 0.827, 0.926, 0.979, 0.897, 0.695, 0.788, 0.816, 0.943, 0.85, 0.761, 0.997, 0.992, 0.963, 0.897, 0.807, 0.781, 0.949, 0.547, 0.964, 0.508, 0.576, 0.64, 0.514, 0.556, 0.777, 0.674, 0.805, 0.727, 0.845, 0.787, 0.889, 0.124, 0.171, 0.47, 0.908, 0.057, 0.106, 0.872, 0.786, 0.86, 0.698, 0.387, 0.985, 0.821, 0.638, 0.545, 0.962, 0.897, 0.926, 0.795, 0.768, 0.692, 0.7, 0.901, 0.275, 0.662, 0.339, 0.946, 0.691, 0.969, 0.917, 0.383, 0.514, 0.315, 0.71, 0.917, 0.882, 0.949, 0.49, 0.802, 0.841, 0.857, 0.919, 0.704, 0.791, 0.774, 0.974, 0.934, 0.583, 0.694, 0.886, 0.921, 0.478, 0.762, 0.999, 0.063, 0.971, 0.938, 0.865, 0.782, 0.413, 0.963, 0.773, 0.909, 0.647, 0.954, 0.782, 0.794, 0.943, 0.395, 0.778, 0.438, 0.874, 0.91, 0.403, 0.882, 0.997, 0.81, 0.6, 0.836, 0.539, 0.498, 0.832, 0.766, 0.727, 0.984, 0.992, 0.425, 0.894, 0.924, 0.997, 0.905, 0.943, 0.668, 0.601, 0.652, 0.984, 0.902, 0.697, 0.497]
global origin = 1
global destination = 50