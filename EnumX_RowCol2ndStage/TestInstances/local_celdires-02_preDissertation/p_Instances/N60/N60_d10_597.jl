global arcs = [1 20; 1 41; 2 20; 2 46; 2 49; 2 60; 3 13; 3 15; 3 45; 3 47; 3 53; 4 6; 4 10; 4 44; 4 47; 5 13; 5 17; 5 22; 5 26; 5 33; 5 46; 5 54; 5 59; 6 8; 6 17; 6 21; 6 23; 6 46; 7 4; 7 10; 7 14; 7 18; 7 23; 7 35; 7 40; 7 43; 7 52; 8 6; 8 10; 8 14; 8 23; 8 29; 8 36; 8 38; 8 52; 8 56; 9 29; 9 37; 9 39; 9 52; 10 30; 10 40; 11 39; 11 42; 12 5; 12 9; 12 19; 12 29; 12 37; 13 6; 13 30; 13 53; 14 19; 14 27; 14 40; 15 11; 15 20; 15 33; 15 39; 15 40; 15 46; 15 57; 15 58; 16 7; 16 13; 16 39; 16 41; 16 53; 16 56; 16 57; 17 45; 17 46; 18 19; 18 24; 18 37; 18 45; 18 53; 18 58; 19 11; 19 12; 19 16; 19 37; 19 41; 19 51; 19 56; 20 4; 20 7; 20 11; 20 17; 20 30; 20 42; 20 44; 20 46; 20 55; 21 3; 21 15; 21 28; 21 39; 21 54; 22 6; 22 29; 22 31; 22 34; 22 39; 23 11; 23 19; 23 32; 23 34; 23 49; 23 55; 23 57; 24 11; 24 29; 24 33; 24 40; 24 48; 25 2; 25 11; 25 13; 25 20; 25 26; 26 19; 26 29; 26 32; 26 36; 26 50; 26 53; 26 54; 26 55; 27 6; 27 25; 27 29; 27 34; 27 50; 28 2; 28 32; 28 58; 29 11; 29 32; 30 10; 30 12; 30 23; 30 31; 30 59; 31 5; 31 19; 31 40; 31 43; 31 59; 32 18; 32 47; 32 53; 33 3; 33 4; 33 8; 33 9; 33 14; 33 15; 33 26; 33 30; 33 31; 33 39; 33 43; 33 46; 33 58; 34 8; 34 16; 34 21; 34 22; 34 27; 34 48; 34 52; 34 60; 35 11; 35 13; 35 15; 35 29; 35 37; 35 47; 35 59; 36 10; 36 24; 36 33; 36 41; 37 11; 37 12; 37 15; 37 21; 37 27; 37 29; 37 42; 37 45; 38 9; 38 13; 38 30; 38 42; 38 53; 39 9; 39 19; 39 28; 39 29; 39 50; 39 58; 40 6; 40 30; 40 32; 40 48; 40 60; 41 21; 41 22; 41 25; 41 33; 41 47; 41 48; 41 51; 41 54; 42 23; 42 24; 42 30; 43 12; 43 21; 43 24; 43 28; 43 40; 43 48; 43 53; 44 2; 44 26; 44 32; 44 52; 44 56; 45 11; 45 28; 45 31; 45 34; 45 37; 45 47; 45 49; 46 4; 46 8; 46 19; 46 24; 46 54; 46 60; 47 10; 47 43; 48 2; 48 24; 48 58; 49 4; 49 19; 49 20; 49 26; 49 47; 50 10; 50 17; 50 28; 50 30; 50 45; 50 46; 50 47; 50 56; 51 12; 51 50; 51 60; 52 22; 52 51; 52 53; 53 4; 53 12; 53 32; 53 36; 53 54; 53 57; 53 59; 53 60; 54 9; 54 21; 54 25; 54 42; 54 44; 54 50; 54 55; 54 59; 55 31; 55 44; 56 13; 56 21; 56 39; 56 41; 56 47; 57 4; 57 27; 57 29; 57 40; 57 42; 58 17; 58 28; 58 30; 58 39; 58 41; 58 56; 59 7; 59 12; 59 18; 59 46; 59 49; 59 51; 59 52]
global d_x = [9.0, 8.0, 3.0, 3.0, 8.0, 4.0, 8.0, 6.0, 2.0, 6.0, 5.0, 6.0, 2.0, 4.0, 1.0, 9.0, 4.0, 10.0, 9.0, 9.0, 8.0, 5.0, 5.0, 6.0, 8.0, 9.0, 3.0, 1.0, 6.0, 9.0, 9.0, 7.0, 9.0, 6.0, 5.0, 3.0, 2.0, 2.0, 7.0, 3.0, 6.0, 10.0, 6.0, 10.0, 2.0, 2.0, 3.0, 10.0, 2.0, 10.0, 4.0, 3.0, 1.0, 3.0, 9.0, 8.0, 4.0, 10.0, 10.0, 1.0, 8.0, 1.0, 1.0, 8.0, 10.0, 1.0, 5.0, 9.0, 9.0, 10.0, 5.0, 8.0, 1.0, 2.0, 4.0, 3.0, 4.0, 2.0, 10.0, 5.0, 3.0, 8.0, 4.0, 1.0, 10.0, 9.0, 2.0, 10.0, 5.0, 7.0, 1.0, 5.0, 10.0, 3.0, 6.0, 5.0, 5.0, 6.0, 10.0, 9.0, 10.0, 8.0, 1.0, 4.0, 8.0, 2.0, 7.0, 8.0, 3.0, 9.0, 7.0, 1.0, 10.0, 1.0, 9.0, 3.0, 7.0, 8.0, 1.0, 5.0, 3.0, 2.0, 3.0, 4.0, 10.0, 7.0, 7.0, 2.0, 2.0, 1.0, 5.0, 6.0, 10.0, 4.0, 9.0, 3.0, 3.0, 4.0, 7.0, 1.0, 8.0, 1.0, 7.0, 5.0, 6.0, 8.0, 9.0, 9.0, 5.0, 10.0, 10.0, 7.0, 5.0, 2.0, 1.0, 6.0, 6.0, 6.0, 2.0, 2.0, 8.0, 8.0, 5.0, 5.0, 10.0, 5.0, 3.0, 9.0, 8.0, 6.0, 3.0, 10.0, 4.0, 7.0, 9.0, 4.0, 2.0, 5.0, 5.0, 3.0, 1.0, 3.0, 3.0, 1.0, 9.0, 9.0, 10.0, 10.0, 10.0, 1.0, 4.0, 8.0, 7.0, 9.0, 4.0, 8.0, 6.0, 2.0, 1.0, 7.0, 4.0, 4.0, 7.0, 2.0, 10.0, 2.0, 8.0, 9.0, 2.0, 7.0, 7.0, 7.0, 10.0, 4.0, 2.0, 9.0, 10.0, 5.0, 10.0, 4.0, 5.0, 10.0, 1.0, 3.0, 2.0, 1.0, 1.0, 2.0, 7.0, 2.0, 7.0, 4.0, 2.0, 1.0, 5.0, 2.0, 3.0, 7.0, 4.0, 8.0, 10.0, 7.0, 10.0, 5.0, 7.0, 7.0, 6.0, 6.0, 3.0, 6.0, 6.0, 1.0, 7.0, 10.0, 4.0, 7.0, 8.0, 3.0, 6.0, 9.0, 2.0, 4.0, 2.0, 1.0, 6.0, 5.0, 2.0, 10.0, 3.0, 1.0, 9.0, 1.0, 7.0, 6.0, 2.0, 8.0, 5.0, 9.0, 9.0, 1.0, 7.0, 4.0, 9.0, 2.0, 2.0, 6.0, 7.0, 8.0, 6.0, 6.0, 4.0, 8.0, 9.0, 9.0, 5.0, 2.0, 1.0, 2.0, 10.0, 10.0, 9.0, 1.0, 9.0, 1.0, 5.0, 6.0, 3.0, 5.0, 3.0, 4.0, 4.0, 8.0, 10.0, 8.0, 3.0, 7.0, 1.0, 4.0, 8.0]
global b_x = 5
global d_y = [7.0, 8.0, 1.0, 2.0, 3.0, 6.0, 4.0, 9.0, 9.0, 10.0, 10.0, 2.0, 2.0, 4.0, 1.0, 4.0, 6.0, 7.0, 7.0, 9.0, 10.0, 2.0, 2.0, 4.0, 8.0, 9.0, 8.0, 4.0, 5.0, 8.0, 2.0, 1.0, 7.0, 1.0, 5.0, 3.0, 6.0, 8.0, 6.0, 7.0, 5.0, 6.0, 5.0, 5.0, 5.0, 3.0, 4.0, 10.0, 5.0, 8.0, 7.0, 6.0, 6.0, 2.0, 5.0, 8.0, 9.0, 3.0, 3.0, 10.0, 8.0, 1.0, 3.0, 10.0, 10.0, 9.0, 10.0, 1.0, 9.0, 6.0, 9.0, 8.0, 10.0, 9.0, 8.0, 7.0, 8.0, 10.0, 4.0, 5.0, 10.0, 8.0, 6.0, 8.0, 5.0, 9.0, 5.0, 1.0, 1.0, 6.0, 1.0, 1.0, 8.0, 6.0, 6.0, 8.0, 3.0, 2.0, 3.0, 7.0, 5.0, 5.0, 10.0, 4.0, 4.0, 7.0, 10.0, 5.0, 7.0, 10.0, 10.0, 6.0, 8.0, 9.0, 10.0, 2.0, 4.0, 3.0, 7.0, 2.0, 8.0, 5.0, 4.0, 6.0, 4.0, 3.0, 6.0, 1.0, 3.0, 2.0, 6.0, 10.0, 10.0, 8.0, 4.0, 1.0, 3.0, 1.0, 9.0, 10.0, 9.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 5.0, 8.0, 3.0, 8.0, 9.0, 8.0, 2.0, 6.0, 2.0, 4.0, 4.0, 9.0, 6.0, 10.0, 9.0, 2.0, 9.0, 5.0, 8.0, 9.0, 1.0, 1.0, 9.0, 9.0, 8.0, 4.0, 10.0, 4.0, 3.0, 8.0, 5.0, 5.0, 5.0, 1.0, 6.0, 1.0, 8.0, 3.0, 7.0, 10.0, 9.0, 8.0, 4.0, 2.0, 1.0, 10.0, 3.0, 4.0, 8.0, 6.0, 6.0, 5.0, 8.0, 2.0, 6.0, 3.0, 3.0, 2.0, 9.0, 5.0, 10.0, 8.0, 6.0, 5.0, 8.0, 5.0, 6.0, 10.0, 4.0, 4.0, 8.0, 10.0, 6.0, 3.0, 1.0, 10.0, 8.0, 8.0, 9.0, 10.0, 5.0, 3.0, 8.0, 9.0, 1.0, 8.0, 1.0, 1.0, 6.0, 5.0, 2.0, 3.0, 5.0, 10.0, 4.0, 1.0, 4.0, 10.0, 7.0, 4.0, 8.0, 1.0, 6.0, 4.0, 6.0, 6.0, 8.0, 1.0, 3.0, 3.0, 10.0, 7.0, 6.0, 6.0, 9.0, 9.0, 10.0, 6.0, 5.0, 8.0, 4.0, 2.0, 7.0, 4.0, 10.0, 8.0, 4.0, 1.0, 2.0, 1.0, 8.0, 9.0, 6.0, 7.0, 6.0, 6.0, 1.0, 7.0, 2.0, 3.0, 1.0, 8.0, 8.0, 3.0, 9.0, 8.0, 7.0, 8.0, 9.0, 8.0, 7.0, 3.0, 9.0, 3.0, 6.0, 8.0, 4.0, 6.0, 3.0, 3.0, 2.0, 9.0, 10.0, 4.0, 7.0, 5.0, 8.0, 5.0, 10.0, 5.0, 5.0, 4.0]
global b_y = 10
global p = [0.143, 0.784, 0.798, 0.512, 0.663, 0.988, 0.881, 0.39, 0.219, 0.805, 0.652, 0.951, 0.614, 0.016, 0.405, 0.979, 0.995, 0.091, 0.825, 0.491, 0.489, 0.971, 0.294, 0.004, 0.176, 0.045, 0.125, 0.678, 0.043, 0.486, 0.56, 0.194, 0.973, 0.785, 0.511, 0.454, 0.385, 0.699, 0.601, 0.617, 0.485, 0.578, 0.377, 0.566, 0.69, 0.712, 0.378, 0.645, 0.393, 0.09, 0.139, 0.59, 0.495, 0.657, 0.36, 0.42, 0.599, 0.846, 0.394, 0.425, 0.27, 0.978, 0.947, 0.065, 0.778, 0.789, 0.397, 0.166, 0.841, 0.314, 0.105, 0.313, 0.042, 0.431, 0.737, 0.467, 0.007, 0.602, 0.562, 0.271, 0.564, 0.447, 0.26, 0.428, 0.232, 0.07, 0.098, 0.093, 0.676, 0.376, 0.732, 0.487, 0.612, 0.011, 0.184, 0.263, 0.995, 0.988, 0.072, 0.696, 0.689, 0.833, 0.206, 0.379, 0.397, 0.419, 0.691, 0.29, 0.795, 0.703, 0.103, 0.319, 0.78, 0.32, 0.057, 0.641, 0.281, 0.866, 0.919, 0.533, 0.964, 0.745, 0.54, 0.432, 0.143, 0.053, 0.783, 0.068, 0.753, 0.358, 0.162, 0.905, 0.571, 0.28, 0.207, 0.4, 0.199, 0.43, 0.046, 0.805, 0.883, 0.923, 0.481, 0.172, 0.83, 0.286, 0.257, 0.019, 0.247, 0.614, 0.048, 0.541, 0.408, 0.866, 0.752, 0.434, 0.348, 0.899, 0.379, 0.948, 0.846, 0.068, 0.639, 0.804, 0.955, 0.942, 0.008, 0.925, 0.407, 0.459, 0.813, 0.948, 0.735, 0.552, 0.134, 0.926, 0.923, 0.476, 0.087, 0.092, 0.228, 0.397, 0.343, 0.573, 0.319, 0.594, 0.739, 0.915, 0.134, 0.941, 0.541, 0.403, 0.494, 0.603, 0.488, 0.997, 0.483, 0.403, 0.331, 0.535, 0.045, 0.408, 0.002, 0.274, 0.944, 0.883, 0.555, 0.076, 0.73, 0.633, 0.517, 0.249, 0.17, 0.348, 0.586, 0.932, 0.525, 0.398, 0.93, 0.26, 0.42, 0.265, 0.664, 0.014, 0.525, 0.521, 0.127, 0.686, 0.473, 0.438, 0.518, 0.56, 0.728, 0.75, 0.223, 0.142, 0.008, 0.035, 0.013, 0.501, 0.145, 0.121, 0.924, 0.554, 0.446, 0.478, 0.304, 0.338, 0.652, 0.294, 0.495, 0.091, 0.198, 0.59, 0.799, 0.386, 0.243, 0.796, 0.461, 0.751, 0.559, 0.099, 0.028, 0.545, 0.972, 0.191, 0.677, 0.451, 0.932, 0.246, 0.56, 0.644, 0.835, 0.649, 0.778, 0.794, 0.547, 0.584, 0.954, 0.46, 0.818, 0.376, 0.865, 0.242, 0.896, 0.44, 0.784, 0.64, 0.697, 0.866, 0.652, 0.991, 0.585, 0.467, 0.095, 0.686, 0.972, 0.519, 0.431, 0.596, 0.682, 0.106, 0.069, 0.833, 0.763, 0.406, 0.485, 0.865, 0.235, 0.512, 0.693, 0.196, 0.445, 0.462, 0.586, 0.884, 0.147, 0.071, 0.036]
global q = [0.985, 0.803, 0.803, 0.9, 0.88, 0.992, 0.912, 0.672, 0.959, 0.995, 0.697, 0.995, 0.97, 0.222, 0.616, 0.999, 0.996, 0.468, 0.832, 0.831, 0.982, 0.998, 0.986, 0.481, 0.501, 0.473, 0.577, 0.845, 0.85, 0.908, 0.642, 0.838, 0.976, 0.935, 0.542, 0.824, 0.432, 0.946, 0.821, 0.939, 0.873, 0.978, 0.623, 0.626, 0.79, 0.718, 0.608, 0.824, 0.522, 0.188, 0.692, 0.69, 0.878, 0.933, 0.729, 0.537, 0.673, 0.927, 0.538, 0.598, 0.737, 0.992, 0.957, 0.221, 0.8, 0.886, 0.876, 0.937, 0.94, 0.578, 0.143, 0.838, 0.969, 0.478, 0.767, 0.822, 0.188, 0.94, 0.976, 0.486, 0.876, 0.877, 0.965, 0.626, 0.364, 0.102, 0.615, 0.585, 0.693, 0.899, 0.998, 0.565, 0.851, 0.588, 0.703, 0.549, 0.999, 0.991, 0.48, 0.723, 0.948, 0.916, 0.528, 0.746, 0.576, 0.533, 0.888, 0.685, 0.818, 0.922, 0.719, 0.345, 0.794, 0.963, 0.957, 0.671, 0.433, 0.945, 0.996, 0.924, 0.994, 0.792, 0.957, 0.529, 0.924, 0.368, 0.883, 0.522, 0.831, 0.938, 0.372, 0.974, 0.837, 0.428, 0.436, 0.422, 0.215, 0.98, 0.939, 0.873, 0.925, 0.965, 0.768, 0.986, 0.854, 0.372, 0.555, 0.415, 0.494, 0.927, 0.121, 0.825, 0.807, 0.945, 0.757, 0.58, 0.963, 0.981, 0.833, 0.964, 0.964, 0.916, 0.834, 0.865, 0.966, 0.996, 0.967, 0.976, 0.877, 0.649, 0.999, 0.991, 0.946, 0.657, 0.165, 0.953, 0.999, 0.908, 0.745, 0.314, 0.808, 0.546, 0.622, 0.637, 0.623, 0.862, 0.752, 0.919, 0.642, 0.998, 0.801, 0.948, 0.891, 0.771, 0.825, 0.998, 0.764, 0.881, 0.773, 0.783, 0.287, 0.683, 0.979, 0.757, 0.976, 0.896, 0.745, 0.526, 0.821, 0.673, 0.63, 0.909, 0.231, 0.529, 0.74, 0.942, 0.829, 0.701, 0.989, 0.605, 0.574, 0.9, 0.729, 0.413, 0.808, 0.747, 0.999, 0.691, 0.771, 0.649, 0.848, 0.827, 0.925, 0.851, 0.889, 0.669, 0.435, 0.619, 0.095, 0.556, 0.696, 0.9, 0.994, 0.792, 0.689, 0.771, 0.394, 0.575, 0.956, 0.693, 0.748, 0.49, 0.797, 0.818, 0.987, 0.659, 0.529, 0.846, 0.773, 0.992, 0.647, 0.542, 0.66, 0.816, 0.996, 0.512, 0.98, 0.985, 0.943, 0.573, 0.956, 0.898, 0.871, 0.956, 0.851, 0.796, 0.676, 0.628, 0.968, 0.73, 0.839, 0.6, 0.938, 0.544, 0.91, 0.908, 0.806, 0.696, 0.727, 0.994, 0.845, 0.994, 0.893, 0.553, 0.646, 0.991, 0.995, 0.902, 0.76, 0.905, 0.756, 0.222, 0.896, 0.839, 0.945, 0.818, 0.503, 0.903, 0.957, 0.874, 0.852, 0.494, 0.56, 0.901, 0.693, 0.961, 0.643, 0.315, 0.149]
global origin = 1
global destination = 60