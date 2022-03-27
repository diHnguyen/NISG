global arcs = [1 3; 1 6; 1 10; 1 26; 1 31; 1 38; 2 15; 2 23; 2 32; 2 55; 2 56; 2 58; 3 8; 3 41; 3 50; 3 59; 4 9; 4 38; 4 50; 5 10; 5 27; 5 37; 6 21; 6 24; 6 27; 6 40; 6 49; 6 58; 7 4; 7 11; 7 25; 7 29; 7 45; 7 46; 7 52; 8 14; 8 22; 8 55; 9 11; 9 27; 9 52; 9 56; 10 7; 10 47; 11 5; 11 35; 11 41; 11 56; 11 57; 12 4; 12 14; 12 27; 12 32; 12 38; 12 41; 12 56; 13 9; 13 19; 13 22; 13 28; 13 38; 13 39; 14 12; 14 30; 14 47; 14 50; 14 54; 14 58; 15 4; 15 8; 15 9; 15 44; 15 46; 15 53; 16 2; 16 10; 16 14; 16 28; 16 42; 16 49; 16 53; 17 9; 17 16; 17 30; 17 34; 17 44; 18 31; 18 36; 18 48; 18 57; 19 4; 19 14; 19 17; 19 37; 19 38; 19 42; 19 54; 20 3; 20 7; 20 30; 20 34; 20 42; 20 55; 21 26; 21 46; 21 48; 22 7; 22 18; 22 46; 22 54; 22 56; 22 58; 22 59; 23 4; 23 24; 23 27; 23 48; 24 10; 24 26; 24 36; 24 57; 24 60; 25 2; 25 13; 25 22; 25 31; 25 41; 26 12; 26 24; 26 45; 26 58; 26 59; 27 7; 27 8; 27 22; 27 23; 27 37; 27 41; 28 20; 28 30; 28 49; 28 53; 28 58; 28 59; 29 3; 29 7; 29 31; 29 44; 29 48; 30 6; 30 11; 30 16; 30 18; 30 33; 30 47; 30 53; 31 4; 31 15; 31 20; 31 24; 31 26; 31 41; 31 42; 31 47; 32 5; 32 17; 32 41; 32 44; 32 59; 33 17; 33 24; 33 44; 33 48; 33 53; 34 2; 34 3; 34 14; 34 15; 34 38; 34 43; 34 46; 34 47; 34 52; 35 5; 35 28; 35 41; 35 47; 35 51; 35 58; 36 11; 36 12; 36 21; 36 27; 36 35; 36 37; 36 54; 36 58; 37 7; 37 14; 37 18; 37 35; 37 43; 37 59; 38 5; 38 6; 38 24; 38 47; 39 6; 39 8; 39 40; 39 43; 39 50; 39 56; 40 4; 40 44; 40 47; 40 49; 40 58; 40 60; 41 22; 41 29; 42 16; 42 29; 42 44; 42 49; 42 54; 43 2; 43 15; 43 21; 43 26; 43 30; 43 32; 43 40; 43 49; 43 52; 43 56; 43 58; 43 59; 44 11; 44 12; 44 52; 44 55; 44 60; 45 15; 45 32; 45 37; 45 39; 45 41; 46 6; 46 12; 46 14; 46 26; 46 36; 46 43; 46 45; 46 47; 47 14; 47 17; 47 36; 47 38; 47 46; 47 56; 48 6; 48 37; 48 38; 49 8; 49 16; 49 17; 49 20; 49 26; 49 27; 49 30; 49 52; 49 53; 50 12; 50 13; 50 28; 50 29; 50 34; 50 38; 50 45; 51 26; 51 31; 51 42; 51 44; 51 48; 52 5; 52 6; 52 15; 52 30; 52 34; 52 36; 52 37; 52 46; 52 51; 52 60; 53 36; 53 46; 53 60; 54 24; 54 33; 54 36; 54 47; 55 3; 55 40; 55 45; 55 48; 56 14; 56 18; 56 49; 56 51; 56 52; 56 54; 57 3; 57 6; 57 12; 57 39; 57 43; 57 44; 58 15; 58 19; 58 22; 58 32; 58 37; 58 40; 58 43; 58 44; 58 56; 59 4; 59 7; 59 9; 59 41; 59 48; 59 60]
global d_x = [4.0, 7.0, 1.0, 9.0, 6.0, 9.0, 10.0, 5.0, 1.0, 7.0, 2.0, 6.0, 6.0, 5.0, 6.0, 3.0, 8.0, 6.0, 7.0, 4.0, 6.0, 8.0, 9.0, 10.0, 6.0, 2.0, 2.0, 10.0, 7.0, 7.0, 5.0, 2.0, 1.0, 1.0, 5.0, 3.0, 8.0, 8.0, 6.0, 5.0, 6.0, 5.0, 9.0, 9.0, 9.0, 7.0, 3.0, 10.0, 5.0, 5.0, 6.0, 6.0, 10.0, 7.0, 7.0, 6.0, 7.0, 5.0, 3.0, 1.0, 1.0, 10.0, 3.0, 2.0, 2.0, 10.0, 9.0, 6.0, 6.0, 2.0, 7.0, 8.0, 7.0, 5.0, 6.0, 3.0, 3.0, 1.0, 2.0, 8.0, 3.0, 5.0, 9.0, 1.0, 7.0, 4.0, 8.0, 3.0, 1.0, 4.0, 3.0, 4.0, 4.0, 9.0, 10.0, 9.0, 10.0, 4.0, 7.0, 5.0, 7.0, 8.0, 8.0, 3.0, 5.0, 5.0, 6.0, 7.0, 3.0, 10.0, 6.0, 7.0, 5.0, 2.0, 10.0, 7.0, 4.0, 6.0, 4.0, 1.0, 7.0, 10.0, 4.0, 2.0, 2.0, 4.0, 5.0, 2.0, 5.0, 7.0, 5.0, 3.0, 6.0, 9.0, 6.0, 5.0, 5.0, 8.0, 10.0, 9.0, 3.0, 4.0, 1.0, 10.0, 1.0, 6.0, 6.0, 10.0, 7.0, 3.0, 8.0, 10.0, 10.0, 7.0, 10.0, 8.0, 1.0, 4.0, 4.0, 5.0, 2.0, 2.0, 9.0, 9.0, 3.0, 6.0, 6.0, 1.0, 9.0, 8.0, 5.0, 7.0, 6.0, 3.0, 9.0, 1.0, 8.0, 8.0, 2.0, 10.0, 1.0, 4.0, 7.0, 3.0, 1.0, 4.0, 8.0, 2.0, 1.0, 8.0, 8.0, 4.0, 7.0, 5.0, 1.0, 9.0, 2.0, 5.0, 7.0, 5.0, 6.0, 9.0, 4.0, 9.0, 3.0, 6.0, 2.0, 4.0, 2.0, 6.0, 5.0, 10.0, 4.0, 1.0, 3.0, 8.0, 9.0, 10.0, 6.0, 6.0, 3.0, 2.0, 5.0, 10.0, 4.0, 6.0, 9.0, 10.0, 1.0, 10.0, 4.0, 6.0, 8.0, 8.0, 8.0, 4.0, 8.0, 9.0, 2.0, 2.0, 9.0, 6.0, 7.0, 3.0, 5.0, 7.0, 5.0, 5.0, 9.0, 6.0, 3.0, 5.0, 3.0, 9.0, 3.0, 3.0, 10.0, 7.0, 1.0, 7.0, 1.0, 8.0, 4.0, 2.0, 8.0, 3.0, 10.0, 1.0, 8.0, 6.0, 4.0, 1.0, 7.0, 3.0, 9.0, 9.0, 8.0, 2.0, 10.0, 5.0, 7.0, 8.0, 7.0, 6.0, 9.0, 10.0, 3.0, 5.0, 7.0, 9.0, 3.0, 3.0, 6.0, 7.0, 2.0, 6.0, 9.0, 9.0, 10.0, 7.0, 4.0, 4.0, 10.0, 8.0, 7.0, 7.0, 7.0, 10.0, 9.0, 4.0, 8.0, 1.0, 9.0, 10.0, 1.0, 5.0, 8.0, 3.0, 5.0, 6.0, 2.0, 2.0, 6.0, 8.0, 4.0, 3.0, 2.0, 4.0, 10.0, 4.0, 4.0, 8.0, 5.0, 1.0]
global b_x = 5
global d_y = [8.0, 5.0, 3.0, 4.0, 7.0, 2.0, 9.0, 8.0, 7.0, 2.0, 5.0, 2.0, 7.0, 5.0, 8.0, 2.0, 2.0, 10.0, 7.0, 9.0, 4.0, 5.0, 9.0, 6.0, 4.0, 7.0, 1.0, 1.0, 5.0, 3.0, 1.0, 5.0, 3.0, 5.0, 10.0, 9.0, 5.0, 1.0, 9.0, 1.0, 3.0, 9.0, 10.0, 8.0, 5.0, 10.0, 7.0, 6.0, 5.0, 8.0, 5.0, 2.0, 2.0, 2.0, 5.0, 1.0, 7.0, 4.0, 6.0, 3.0, 8.0, 5.0, 4.0, 6.0, 2.0, 2.0, 5.0, 5.0, 8.0, 1.0, 2.0, 3.0, 3.0, 8.0, 5.0, 6.0, 9.0, 9.0, 5.0, 9.0, 8.0, 5.0, 8.0, 7.0, 10.0, 6.0, 8.0, 10.0, 9.0, 5.0, 5.0, 6.0, 7.0, 7.0, 10.0, 7.0, 3.0, 7.0, 6.0, 1.0, 5.0, 6.0, 3.0, 8.0, 1.0, 3.0, 2.0, 7.0, 9.0, 7.0, 1.0, 5.0, 9.0, 8.0, 9.0, 2.0, 1.0, 1.0, 4.0, 6.0, 3.0, 8.0, 2.0, 5.0, 9.0, 5.0, 3.0, 5.0, 9.0, 8.0, 3.0, 4.0, 1.0, 7.0, 8.0, 5.0, 2.0, 2.0, 7.0, 6.0, 10.0, 8.0, 8.0, 7.0, 3.0, 1.0, 1.0, 6.0, 10.0, 3.0, 4.0, 6.0, 10.0, 4.0, 5.0, 2.0, 4.0, 1.0, 9.0, 5.0, 9.0, 5.0, 6.0, 9.0, 9.0, 10.0, 2.0, 9.0, 6.0, 6.0, 3.0, 10.0, 4.0, 8.0, 6.0, 9.0, 10.0, 6.0, 4.0, 9.0, 7.0, 4.0, 6.0, 5.0, 8.0, 6.0, 10.0, 9.0, 1.0, 3.0, 8.0, 8.0, 10.0, 4.0, 6.0, 5.0, 2.0, 7.0, 1.0, 10.0, 3.0, 5.0, 3.0, 9.0, 5.0, 5.0, 3.0, 10.0, 6.0, 2.0, 9.0, 8.0, 5.0, 7.0, 10.0, 2.0, 5.0, 1.0, 8.0, 7.0, 8.0, 5.0, 5.0, 2.0, 5.0, 9.0, 1.0, 6.0, 8.0, 5.0, 8.0, 6.0, 1.0, 2.0, 2.0, 1.0, 1.0, 7.0, 9.0, 8.0, 8.0, 10.0, 4.0, 10.0, 4.0, 9.0, 1.0, 6.0, 7.0, 7.0, 8.0, 2.0, 8.0, 8.0, 6.0, 9.0, 9.0, 9.0, 7.0, 4.0, 10.0, 3.0, 5.0, 8.0, 2.0, 10.0, 8.0, 3.0, 8.0, 7.0, 10.0, 1.0, 8.0, 4.0, 5.0, 6.0, 4.0, 10.0, 1.0, 1.0, 4.0, 2.0, 9.0, 9.0, 2.0, 5.0, 2.0, 7.0, 5.0, 5.0, 3.0, 6.0, 3.0, 6.0, 7.0, 4.0, 4.0, 1.0, 9.0, 2.0, 5.0, 5.0, 6.0, 2.0, 6.0, 2.0, 8.0, 3.0, 1.0, 4.0, 7.0, 7.0, 7.0, 6.0, 1.0, 2.0, 10.0, 6.0, 1.0, 9.0, 8.0, 2.0, 9.0, 5.0, 1.0, 10.0, 7.0, 8.0, 5.0, 9.0, 10.0, 2.0, 7.0, 6.0]
global b_y = 10
global p = [0.095, 0.731, 0.032, 0.948, 0.963, 0.755, 0.279, 0.831, 0.726, 0.183, 0.292, 0.232, 0.497, 0.964, 0.212, 0.277, 0.318, 0.615, 0.054, 0.286, 0.182, 0.6, 0.817, 0.583, 0.156, 0.455, 0.327, 0.318, 0.395, 0.166, 0.614, 0.872, 0.072, 0.25, 0.674, 0.943, 0.995, 0.696, 0.712, 0.955, 0.785, 0.882, 0.441, 0.552, 0.737, 0.314, 0.175, 0.172, 0.631, 0.455, 0.852, 0.244, 0.371, 0.684, 0.507, 0.041, 0.035, 0.059, 0.847, 0.832, 0.97, 0.058, 0.633, 0.108, 0.233, 0.608, 0.519, 0.801, 0.769, 0.294, 0.316, 0.122, 0.861, 0.718, 0.828, 0.132, 0.38, 0.175, 0.765, 0.32, 0.69, 0.318, 0.397, 0.503, 0.847, 0.036, 0.15, 0.141, 0.938, 0.932, 0.489, 0.916, 0.091, 0.107, 0.626, 0.924, 0.035, 0.78, 0.72, 0.628, 0.73, 0.809, 0.472, 0.021, 0.655, 0.872, 0.955, 0.533, 0.756, 0.136, 0.682, 0.327, 0.567, 0.345, 0.029, 0.976, 0.326, 0.576, 0.346, 0.293, 0.466, 0.252, 0.882, 0.592, 0.696, 0.609, 0.782, 0.428, 0.679, 0.172, 0.885, 0.241, 0.049, 0.861, 0.659, 0.642, 0.743, 0.418, 0.211, 0.262, 0.967, 0.222, 0.038, 0.477, 0.793, 0.401, 0.72, 0.748, 0.305, 0.547, 0.452, 0.094, 0.694, 0.221, 0.385, 0.322, 0.848, 0.746, 0.927, 0.651, 0.592, 0.759, 0.8, 0.276, 0.023, 0.778, 0.441, 0.287, 0.925, 0.686, 0.919, 0.855, 0.743, 0.658, 0.514, 0.319, 0.583, 0.728, 0.404, 0.088, 0.818, 0.15, 0.124, 0.685, 0.321, 0.812, 0.372, 0.468, 0.783, 0.887, 0.664, 0.94, 0.296, 0.795, 0.472, 0.05, 0.025, 0.055, 0.038, 0.073, 0.196, 0.937, 0.596, 0.697, 0.744, 0.586, 0.865, 0.013, 0.954, 0.463, 0.519, 0.579, 0.323, 0.49, 0.672, 0.378, 0.79, 0.41, 0.514, 0.262, 0.616, 0.628, 0.816, 0.86, 0.278, 0.586, 0.166, 0.105, 0.892, 0.453, 0.991, 0.449, 0.612, 0.755, 0.132, 0.564, 0.569, 0.811, 0.489, 0.612, 0.238, 0.669, 0.692, 0.35, 0.662, 0.782, 0.774, 0.614, 0.867, 0.322, 0.286, 0.8, 0.5, 0.215, 0.94, 0.592, 0.223, 0.289, 0.97, 0.029, 0.181, 0.063, 0.704, 0.808, 0.776, 0.704, 0.804, 0.359, 0.388, 0.145, 0.537, 0.454, 0.142, 0.9, 0.037, 0.154, 0.145, 0.341, 0.497, 0.365, 0.6, 0.945, 0.965, 0.694, 0.007, 0.437, 0.817, 0.641, 0.343, 0.446, 0.906, 0.474, 0.437, 0.51, 0.722, 0.477, 0.109, 0.763, 0.893, 0.483, 0.848, 0.916, 0.768, 0.993, 0.552, 0.218, 0.743, 0.087, 0.532, 0.508, 0.426, 0.689, 0.463, 0.295, 0.135, 0.694, 0.041, 0.484, 0.549, 0.834, 0.395, 0.05, 0.482, 0.749, 0.89, 0.37, 0.62, 0.333, 0.671, 0.265, 0.22, 0.893, 0.339, 0.892]
global q = [0.505, 0.931, 0.681, 0.974, 0.966, 0.979, 0.37, 0.886, 0.958, 0.913, 0.459, 0.481, 0.995, 0.98, 0.55, 0.483, 0.743, 0.753, 0.964, 0.787, 0.618, 0.638, 0.924, 0.692, 0.265, 0.988, 0.715, 0.805, 0.848, 0.848, 0.777, 0.999, 0.274, 0.361, 0.763, 0.993, 0.998, 0.873, 0.933, 0.981, 0.897, 0.886, 0.888, 0.921, 0.909, 0.611, 0.931, 0.775, 0.706, 0.571, 0.976, 0.726, 0.776, 0.885, 0.75, 0.505, 0.297, 0.345, 0.92, 0.995, 0.99, 0.205, 0.75, 0.921, 0.631, 0.805, 0.955, 0.983, 0.954, 0.881, 0.932, 0.784, 0.893, 0.87, 0.93, 0.39, 0.466, 0.502, 0.91, 0.844, 0.877, 0.436, 0.519, 0.944, 0.864, 0.568, 0.442, 0.914, 0.978, 0.999, 0.708, 0.967, 0.703, 0.787, 0.751, 0.924, 0.934, 0.914, 0.898, 0.744, 0.937, 0.88, 0.643, 0.132, 0.835, 0.947, 0.978, 0.709, 0.885, 0.45, 0.754, 0.668, 0.677, 0.433, 0.665, 0.986, 0.489, 0.958, 0.491, 0.921, 0.877, 0.422, 0.907, 0.738, 0.965, 0.782, 0.866, 0.458, 0.679, 0.495, 0.889, 0.897, 0.06, 0.946, 0.766, 0.869, 0.986, 0.759, 0.884, 0.621, 0.998, 0.286, 0.912, 0.573, 0.963, 0.456, 0.958, 0.992, 0.99, 0.92, 0.666, 0.452, 0.861, 0.904, 0.576, 0.598, 0.85, 0.905, 0.942, 0.657, 0.819, 0.763, 0.891, 0.768, 0.991, 0.854, 0.832, 0.446, 0.988, 0.735, 0.943, 0.934, 0.816, 0.799, 0.724, 0.786, 0.944, 0.821, 0.877, 0.554, 0.827, 0.878, 0.478, 0.786, 0.47, 0.827, 0.679, 0.669, 0.81, 0.899, 0.719, 0.976, 0.325, 0.966, 0.49, 0.4, 0.142, 0.846, 0.847, 0.823, 0.274, 0.966, 0.946, 0.722, 0.839, 0.724, 0.956, 0.504, 0.976, 0.701, 0.797, 0.962, 0.463, 0.933, 0.852, 0.915, 0.803, 0.571, 0.792, 0.561, 0.709, 0.756, 0.969, 0.899, 0.573, 0.973, 0.676, 0.296, 0.927, 0.544, 0.991, 0.701, 0.762, 0.863, 0.891, 0.935, 0.708, 0.946, 0.619, 0.746, 0.41, 0.838, 0.883, 0.907, 0.96, 0.952, 0.935, 0.804, 0.923, 0.847, 0.586, 0.905, 0.542, 0.484, 0.964, 0.704, 0.564, 0.711, 0.994, 0.495, 0.355, 0.804, 0.73, 0.838, 0.847, 0.761, 0.979, 0.938, 0.702, 0.689, 0.839, 0.569, 0.625, 0.911, 0.736, 0.277, 0.954, 0.478, 0.544, 0.901, 0.786, 0.967, 0.98, 0.721, 0.147, 0.545, 0.873, 0.999, 0.933, 0.638, 0.989, 0.519, 0.524, 0.937, 0.885, 0.636, 0.9, 0.795, 0.988, 0.999, 0.88, 0.959, 0.821, 0.999, 0.942, 0.941, 0.749, 0.514, 0.952, 0.731, 0.444, 0.84, 0.835, 0.325, 0.56, 0.915, 0.612, 0.837, 0.992, 0.854, 0.405, 0.531, 0.829, 0.972, 0.995, 0.794, 0.772, 0.452, 0.943, 0.324, 0.345, 0.92, 0.856, 0.975]
global origin = 1
global destination = 60