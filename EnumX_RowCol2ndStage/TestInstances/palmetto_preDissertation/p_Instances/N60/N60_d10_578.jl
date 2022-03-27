global arcs = [1 9; 1 14; 1 19; 1 38; 1 42; 1 54; 1 58; 2 5; 2 25; 2 29; 2 31; 2 32; 2 37; 2 42; 2 43; 2 53; 2 54; 2 58; 3 13; 3 16; 3 32; 3 48; 3 53; 3 55; 3 59; 4 6; 4 34; 4 52; 4 56; 5 13; 5 32; 6 23; 6 29; 6 46; 6 53; 7 3; 7 8; 7 17; 7 23; 7 33; 7 43; 8 9; 8 17; 8 21; 8 22; 8 23; 8 25; 8 27; 8 30; 8 32; 8 44; 8 47; 8 48; 8 51; 8 53; 9 14; 9 22; 9 29; 9 38; 9 51; 10 16; 10 20; 10 34; 11 9; 12 5; 12 11; 12 13; 12 15; 12 21; 12 24; 12 39; 12 40; 12 53; 12 58; 12 60; 13 12; 13 21; 13 24; 13 31; 14 2; 14 9; 14 12; 14 15; 14 16; 14 24; 14 37; 14 39; 14 60; 15 2; 15 3; 15 6; 15 17; 15 60; 16 5; 16 18; 16 31; 16 41; 16 45; 17 4; 17 30; 17 33; 17 47; 17 48; 18 14; 18 30; 18 33; 18 42; 18 48; 18 56; 19 3; 19 10; 19 21; 19 22; 19 23; 19 27; 19 32; 19 49; 19 54; 20 26; 20 30; 20 44; 20 49; 21 5; 21 31; 21 38; 21 41; 21 50; 21 59; 22 11; 22 17; 22 43; 22 47; 23 10; 23 49; 23 59; 23 60; 24 49; 24 53; 25 8; 25 16; 25 17; 25 20; 25 23; 25 24; 25 40; 25 42; 25 49; 26 5; 26 27; 26 29; 26 31; 27 5; 27 60; 28 3; 28 5; 28 12; 28 27; 28 44; 29 6; 29 23; 29 36; 29 37; 29 48; 29 54; 30 2; 30 6; 30 14; 30 16; 30 17; 30 43; 30 60; 31 6; 31 12; 31 13; 31 14; 31 46; 31 59; 32 16; 32 19; 32 31; 32 33; 32 36; 32 51; 32 56; 33 9; 33 36; 33 55; 34 5; 34 9; 34 18; 34 28; 34 30; 35 24; 35 27; 35 45; 35 59; 36 6; 36 17; 36 25; 36 41; 36 48; 36 52; 37 20; 37 32; 37 44; 38 5; 38 11; 38 14; 38 55; 39 15; 39 20; 39 26; 39 32; 39 33; 39 36; 39 42; 40 10; 40 17; 40 24; 40 39; 40 51; 41 20; 41 24; 41 27; 41 33; 41 51; 42 4; 42 10; 42 14; 42 60; 43 7; 43 20; 43 35; 43 40; 43 47; 43 54; 43 55; 44 3; 44 7; 44 22; 44 24; 44 30; 44 32; 44 60; 45 3; 45 12; 45 14; 45 22; 45 34; 45 39; 45 55; 45 59; 46 2; 46 4; 46 5; 46 32; 46 40; 46 54; 47 28; 47 30; 47 33; 47 50; 48 7; 48 29; 48 30; 48 35; 48 45; 48 47; 48 58; 49 3; 49 9; 49 46; 49 48; 49 51; 50 15; 50 20; 50 23; 50 25; 50 28; 50 32; 50 44; 50 49; 51 36; 51 40; 51 55; 51 59; 52 2; 52 4; 52 8; 52 34; 52 42; 53 11; 53 13; 53 17; 53 21; 53 36; 53 44; 53 47; 53 48; 53 52; 53 57; 53 59; 54 30; 54 36; 54 56; 55 8; 55 10; 55 22; 55 31; 55 35; 55 43; 56 8; 56 21; 56 28; 56 31; 56 44; 56 49; 57 10; 57 31; 57 48; 57 56; 58 15; 58 47; 58 48; 59 8; 59 13; 59 15; 59 17; 59 52]
global d_x = [2.0, 7.0, 8.0, 2.0, 9.0, 9.0, 1.0, 7.0, 1.0, 3.0, 6.0, 6.0, 2.0, 2.0, 5.0, 10.0, 9.0, 10.0, 3.0, 6.0, 2.0, 9.0, 7.0, 9.0, 7.0, 5.0, 9.0, 8.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 9.0, 2.0, 6.0, 3.0, 1.0, 10.0, 6.0, 10.0, 1.0, 5.0, 4.0, 5.0, 10.0, 9.0, 3.0, 7.0, 10.0, 5.0, 9.0, 5.0, 10.0, 1.0, 8.0, 2.0, 2.0, 2.0, 10.0, 9.0, 6.0, 3.0, 4.0, 4.0, 7.0, 3.0, 2.0, 1.0, 7.0, 5.0, 4.0, 6.0, 5.0, 7.0, 9.0, 9.0, 10.0, 5.0, 10.0, 8.0, 5.0, 5.0, 6.0, 2.0, 10.0, 8.0, 1.0, 2.0, 1.0, 10.0, 10.0, 4.0, 2.0, 8.0, 7.0, 7.0, 1.0, 4.0, 9.0, 9.0, 4.0, 6.0, 8.0, 4.0, 5.0, 5.0, 7.0, 2.0, 4.0, 10.0, 10.0, 8.0, 7.0, 2.0, 7.0, 9.0, 5.0, 9.0, 5.0, 7.0, 7.0, 6.0, 3.0, 5.0, 6.0, 2.0, 6.0, 2.0, 3.0, 1.0, 5.0, 10.0, 2.0, 10.0, 10.0, 3.0, 10.0, 1.0, 4.0, 1.0, 5.0, 1.0, 10.0, 5.0, 4.0, 5.0, 7.0, 2.0, 9.0, 3.0, 10.0, 6.0, 9.0, 1.0, 1.0, 1.0, 5.0, 8.0, 2.0, 5.0, 8.0, 5.0, 7.0, 5.0, 8.0, 7.0, 3.0, 5.0, 1.0, 10.0, 5.0, 3.0, 3.0, 4.0, 6.0, 2.0, 2.0, 7.0, 1.0, 5.0, 6.0, 7.0, 3.0, 2.0, 2.0, 10.0, 4.0, 3.0, 5.0, 6.0, 10.0, 8.0, 7.0, 5.0, 8.0, 4.0, 9.0, 5.0, 10.0, 9.0, 3.0, 10.0, 8.0, 1.0, 5.0, 1.0, 10.0, 5.0, 2.0, 4.0, 5.0, 3.0, 4.0, 5.0, 1.0, 4.0, 1.0, 9.0, 5.0, 1.0, 6.0, 8.0, 10.0, 3.0, 5.0, 5.0, 1.0, 7.0, 8.0, 3.0, 5.0, 8.0, 3.0, 3.0, 1.0, 8.0, 6.0, 2.0, 8.0, 6.0, 10.0, 10.0, 10.0, 5.0, 7.0, 9.0, 4.0, 5.0, 4.0, 6.0, 8.0, 1.0, 6.0, 5.0, 10.0, 8.0, 9.0, 6.0, 10.0, 8.0, 9.0, 7.0, 10.0, 4.0, 3.0, 10.0, 10.0, 1.0, 9.0, 10.0, 9.0, 10.0, 4.0, 9.0, 5.0, 7.0, 1.0, 2.0, 2.0, 4.0, 3.0, 5.0, 8.0, 2.0, 8.0, 5.0, 3.0, 10.0, 1.0, 7.0, 8.0, 10.0, 9.0, 2.0, 6.0, 9.0, 6.0, 7.0, 3.0, 4.0, 3.0, 3.0, 5.0, 4.0, 4.0, 7.0, 2.0, 8.0, 4.0, 10.0, 1.0, 6.0, 3.0, 3.0, 7.0, 10.0, 10.0, 5.0, 7.0, 7.0, 8.0, 10.0, 5.0, 8.0, 5.0, 6.0, 3.0]
global b_x = 5
global d_y = [1.0, 5.0, 8.0, 10.0, 9.0, 3.0, 10.0, 10.0, 7.0, 10.0, 2.0, 5.0, 3.0, 7.0, 6.0, 6.0, 6.0, 5.0, 9.0, 7.0, 1.0, 10.0, 8.0, 10.0, 7.0, 8.0, 8.0, 10.0, 8.0, 3.0, 8.0, 4.0, 7.0, 10.0, 5.0, 9.0, 7.0, 4.0, 6.0, 9.0, 7.0, 10.0, 9.0, 10.0, 6.0, 7.0, 9.0, 6.0, 4.0, 10.0, 1.0, 2.0, 4.0, 9.0, 1.0, 1.0, 9.0, 9.0, 8.0, 3.0, 9.0, 10.0, 2.0, 6.0, 2.0, 9.0, 8.0, 10.0, 8.0, 10.0, 2.0, 8.0, 9.0, 6.0, 8.0, 8.0, 8.0, 2.0, 7.0, 7.0, 9.0, 6.0, 3.0, 10.0, 9.0, 2.0, 8.0, 5.0, 6.0, 4.0, 8.0, 8.0, 2.0, 6.0, 8.0, 1.0, 9.0, 3.0, 6.0, 9.0, 9.0, 8.0, 3.0, 9.0, 10.0, 4.0, 9.0, 10.0, 9.0, 4.0, 5.0, 6.0, 6.0, 10.0, 1.0, 3.0, 7.0, 6.0, 10.0, 7.0, 3.0, 10.0, 6.0, 2.0, 10.0, 9.0, 2.0, 3.0, 4.0, 10.0, 10.0, 10.0, 1.0, 9.0, 10.0, 7.0, 10.0, 10.0, 1.0, 7.0, 5.0, 3.0, 9.0, 9.0, 4.0, 6.0, 4.0, 4.0, 10.0, 7.0, 7.0, 10.0, 7.0, 1.0, 6.0, 3.0, 7.0, 3.0, 5.0, 2.0, 1.0, 2.0, 4.0, 6.0, 4.0, 10.0, 10.0, 9.0, 5.0, 5.0, 4.0, 9.0, 8.0, 9.0, 5.0, 8.0, 7.0, 3.0, 7.0, 9.0, 7.0, 5.0, 6.0, 1.0, 6.0, 6.0, 3.0, 6.0, 10.0, 4.0, 7.0, 1.0, 1.0, 10.0, 9.0, 3.0, 7.0, 6.0, 8.0, 7.0, 6.0, 8.0, 10.0, 8.0, 1.0, 4.0, 1.0, 9.0, 4.0, 4.0, 10.0, 10.0, 10.0, 4.0, 1.0, 1.0, 8.0, 2.0, 8.0, 4.0, 7.0, 6.0, 10.0, 3.0, 1.0, 9.0, 1.0, 3.0, 3.0, 7.0, 10.0, 9.0, 2.0, 9.0, 8.0, 3.0, 8.0, 2.0, 5.0, 8.0, 1.0, 9.0, 10.0, 4.0, 7.0, 9.0, 9.0, 1.0, 2.0, 2.0, 4.0, 3.0, 3.0, 10.0, 1.0, 6.0, 6.0, 3.0, 9.0, 6.0, 7.0, 7.0, 2.0, 3.0, 1.0, 5.0, 4.0, 3.0, 10.0, 9.0, 7.0, 7.0, 2.0, 5.0, 5.0, 7.0, 1.0, 7.0, 7.0, 9.0, 2.0, 3.0, 5.0, 6.0, 4.0, 9.0, 4.0, 4.0, 1.0, 8.0, 2.0, 7.0, 6.0, 10.0, 8.0, 9.0, 6.0, 7.0, 4.0, 6.0, 3.0, 10.0, 2.0, 3.0, 7.0, 8.0, 3.0, 1.0, 4.0, 6.0, 3.0, 5.0, 8.0, 3.0, 8.0, 10.0, 6.0, 3.0, 3.0, 3.0, 4.0, 3.0, 9.0, 5.0, 5.0, 5.0, 1.0, 5.0, 4.0]
global b_y = 10
global p = [0.98, 0.711, 0.521, 0.575, 0.169, 0.272, 0.412, 0.795, 0.669, 0.837, 0.375, 0.111, 0.911, 0.951, 0.998, 0.236, 0.69, 0.193, 0.209, 0.968, 0.863, 0.559, 0.042, 0.365, 0.975, 0.966, 0.491, 0.606, 0.004, 0.869, 0.885, 0.34, 0.898, 0.623, 0.541, 0.711, 0.071, 0.023, 0.99, 0.722, 0.973, 0.558, 0.87, 0.546, 0.547, 0.019, 0.053, 0.088, 0.325, 0.34, 0.517, 0.765, 0.915, 0.203, 0.722, 0.147, 0.456, 0.349, 0.647, 0.333, 0.79, 0.987, 0.479, 0.279, 0.288, 0.637, 0.756, 0.082, 0.378, 0.878, 0.795, 0.49, 0.92, 0.565, 0.895, 0.577, 0.808, 0.46, 0.338, 0.665, 0.559, 0.474, 0.829, 0.035, 0.736, 0.007, 0.028, 0.952, 0.236, 0.498, 0.927, 0.623, 0.906, 0.654, 0.407, 0.37, 0.197, 0.62, 0.249, 0.365, 0.087, 0.537, 0.138, 0.266, 0.676, 0.161, 0.641, 0.23, 0.508, 0.838, 0.125, 0.604, 0.779, 0.578, 0.709, 0.469, 0.069, 0.165, 0.896, 0.911, 0.83, 0.595, 0.053, 0.982, 0.15, 0.653, 0.336, 0.57, 0.298, 0.67, 0.418, 0.345, 0.342, 0.309, 0.635, 0.178, 0.729, 0.278, 0.917, 0.611, 0.767, 0.882, 0.785, 0.599, 0.007, 0.162, 0.693, 0.555, 0.275, 0.546, 0.17, 0.965, 0.564, 0.5, 0.334, 0.483, 0.848, 0.045, 0.307, 0.92, 0.27, 0.679, 0.168, 0.377, 0.463, 0.007, 0.485, 0.498, 0.374, 0.769, 0.022, 0.589, 0.609, 0.509, 0.341, 0.601, 0.077, 0.741, 0.662, 0.086, 0.864, 0.555, 0.619, 0.379, 0.664, 0.064, 0.475, 0.196, 0.78, 0.251, 0.261, 0.078, 0.35, 0.896, 0.066, 0.88, 0.453, 0.076, 0.204, 0.208, 0.485, 0.657, 0.985, 0.909, 0.574, 0.891, 0.559, 0.501, 0.568, 0.287, 0.607, 0.6, 0.522, 0.665, 0.607, 0.513, 0.226, 0.802, 0.589, 0.535, 0.207, 0.121, 0.466, 0.512, 0.787, 0.786, 0.92, 0.072, 0.707, 0.831, 0.622, 0.826, 0.563, 0.125, 0.433, 0.604, 0.81, 0.008, 0.977, 0.877, 0.601, 0.447, 0.442, 0.805, 0.797, 0.877, 0.027, 0.289, 0.541, 0.641, 0.328, 0.369, 0.421, 0.002, 0.426, 0.5, 0.557, 0.186, 0.743, 0.58, 0.717, 0.172, 0.105, 0.612, 0.407, 0.43, 0.77, 0.744, 0.358, 0.538, 0.639, 0.703, 0.571, 0.672, 0.481, 0.536, 0.376, 0.002, 0.757, 0.52, 0.276, 0.72, 0.184, 0.379, 0.288, 0.658, 0.207, 0.188, 0.036, 0.535, 0.478, 0.097, 0.763, 0.025, 0.715, 0.697, 0.274, 0.652, 0.317, 0.895, 0.491, 0.505, 0.721, 0.993, 0.12, 0.764, 0.83, 0.11, 0.894, 0.121, 0.043, 0.173, 0.525, 0.084, 0.729, 0.919, 0.753, 0.198, 0.025, 0.388, 0.234, 0.521, 0.349, 0.627, 0.631, 0.048, 0.535, 0.197, 0.861]
global q = [0.987, 0.958, 0.878, 0.769, 0.575, 0.86, 0.737, 0.861, 0.803, 0.924, 0.554, 0.545, 0.986, 0.978, 0.998, 0.276, 0.809, 0.851, 0.427, 0.985, 0.985, 0.894, 0.52, 0.832, 0.988, 0.977, 0.681, 0.793, 0.747, 0.992, 0.992, 0.628, 0.935, 0.756, 0.944, 0.759, 0.22, 0.622, 0.994, 0.793, 0.979, 0.704, 0.942, 0.943, 0.951, 0.892, 0.827, 0.667, 0.844, 0.358, 0.889, 0.95, 0.967, 0.334, 0.984, 0.418, 0.483, 0.736, 0.958, 0.42, 0.845, 0.987, 0.857, 0.401, 0.723, 0.719, 0.77, 0.949, 0.381, 0.968, 0.807, 0.806, 0.94, 0.844, 0.984, 0.844, 0.931, 0.718, 0.907, 0.81, 0.677, 0.907, 0.952, 0.101, 0.984, 0.892, 0.086, 0.984, 0.903, 0.778, 0.994, 0.974, 0.986, 0.919, 0.652, 0.847, 0.675, 0.637, 0.385, 0.694, 0.286, 0.962, 0.313, 0.411, 0.736, 0.217, 0.83, 0.734, 0.853, 0.898, 0.886, 0.671, 0.957, 0.845, 0.775, 0.993, 0.987, 0.517, 0.961, 0.973, 0.912, 0.636, 0.846, 0.992, 0.394, 0.719, 0.692, 0.636, 0.437, 0.998, 0.798, 0.718, 0.622, 0.457, 0.735, 0.205, 0.769, 0.907, 0.938, 0.964, 0.885, 0.992, 0.995, 0.993, 0.735, 0.781, 0.753, 0.565, 0.753, 0.631, 0.336, 0.968, 0.966, 0.814, 0.893, 0.583, 0.964, 0.718, 0.37, 0.955, 0.916, 0.714, 0.351, 0.622, 0.826, 0.85, 0.692, 0.607, 0.497, 0.866, 0.479, 0.847, 0.905, 0.574, 0.789, 0.956, 0.587, 0.802, 0.935, 0.262, 0.88, 0.992, 0.735, 0.517, 0.797, 0.241, 0.879, 0.259, 0.838, 0.259, 0.356, 0.3, 0.604, 0.947, 0.13, 0.953, 0.808, 0.234, 0.451, 0.832, 0.867, 0.69, 0.987, 0.998, 0.777, 0.963, 0.583, 0.606, 0.836, 0.672, 0.926, 0.704, 0.677, 0.701, 0.95, 0.838, 0.962, 0.895, 0.994, 0.805, 0.73, 0.19, 0.976, 0.701, 0.845, 0.805, 0.932, 0.941, 0.754, 0.946, 0.695, 0.947, 0.86, 0.202, 0.831, 0.709, 0.937, 0.989, 0.996, 0.892, 0.907, 0.568, 0.479, 0.89, 0.933, 0.921, 0.687, 0.768, 0.879, 0.696, 0.864, 0.927, 0.454, 0.581, 0.763, 0.765, 0.644, 0.465, 0.899, 0.956, 0.844, 0.956, 0.482, 0.864, 0.524, 0.962, 0.836, 0.996, 0.675, 0.657, 0.901, 0.845, 0.699, 0.755, 0.855, 0.545, 0.899, 0.687, 0.99, 0.731, 0.883, 0.882, 0.297, 0.894, 0.837, 0.801, 0.878, 0.656, 0.363, 0.777, 0.515, 0.412, 0.833, 0.48, 0.842, 0.849, 0.39, 0.739, 0.878, 0.899, 0.79, 0.77, 0.891, 0.999, 0.123, 0.979, 0.934, 0.257, 0.924, 0.384, 0.932, 0.702, 0.561, 0.817, 0.96, 0.931, 0.862, 0.817, 0.751, 0.841, 0.59, 0.909, 0.933, 0.699, 0.711, 0.162, 0.951, 0.793, 0.906]
global origin = 1
global destination = 60