global arcs = [1 6; 1 8; 1 12; 1 18; 1 35; 1 38; 2 12; 2 20; 2 22; 2 49; 3 28; 3 30; 3 39; 3 46; 4 22; 4 27; 4 43; 4 49; 5 8; 5 9; 5 19; 5 29; 5 35; 5 43; 5 45; 6 2; 6 5; 6 20; 6 21; 6 32; 6 44; 7 3; 7 6; 7 17; 7 27; 7 28; 7 29; 7 30; 7 35; 7 47; 8 2; 8 5; 8 13; 8 18; 8 44; 9 7; 9 13; 9 21; 9 30; 9 41; 9 44; 9 48; 10 24; 10 39; 10 43; 10 48; 11 35; 11 42; 12 3; 12 6; 12 8; 12 14; 12 28; 12 33; 13 11; 13 18; 13 29; 13 41; 13 46; 14 6; 14 10; 14 25; 14 28; 14 29; 14 37; 14 50; 15 24; 15 33; 15 39; 15 45; 16 21; 16 23; 16 24; 16 32; 17 7; 17 11; 17 21; 17 32; 17 38; 18 7; 18 12; 18 13; 18 16; 18 21; 18 24; 18 25; 18 37; 18 41; 19 11; 19 20; 19 37; 19 48; 20 13; 20 29; 20 49; 21 5; 21 26; 21 31; 21 32; 21 44; 21 47; 21 49; 21 50; 22 7; 22 16; 22 17; 22 34; 22 35; 22 41; 22 48; 23 8; 23 25; 23 40; 24 15; 24 19; 24 23; 24 29; 24 34; 24 47; 24 49; 24 50; 25 22; 25 24; 25 27; 25 35; 25 42; 26 4; 26 8; 26 10; 26 21; 26 24; 26 31; 26 46; 27 6; 27 9; 27 50; 28 13; 28 27; 28 32; 28 36; 28 45; 28 46; 29 8; 29 10; 29 34; 29 44; 29 47; 30 8; 30 28; 30 33; 30 38; 30 43; 31 21; 31 22; 31 34; 32 6; 32 7; 32 19; 32 21; 32 28; 32 31; 32 37; 33 18; 33 38; 33 39; 34 44; 35 32; 36 3; 36 8; 36 9; 36 15; 36 32; 36 33; 36 38; 36 41; 37 6; 37 10; 37 14; 37 28; 37 45; 38 3; 38 4; 38 12; 38 19; 38 30; 38 35; 38 45; 38 48; 39 7; 39 13; 39 31; 39 34; 40 3; 40 9; 40 12; 40 15; 40 18; 40 23; 40 25; 40 46; 41 7; 41 14; 41 40; 41 49; 42 7; 42 21; 42 26; 42 29; 42 31; 42 43; 43 9; 43 10; 43 13; 43 28; 43 29; 43 30; 43 35; 43 46; 43 47; 43 49; 43 50; 44 6; 44 9; 44 17; 44 18; 44 27; 45 11; 45 28; 45 35; 45 47; 45 49; 46 6; 46 11; 46 14; 46 23; 46 30; 46 31; 46 36; 46 49; 47 11; 47 25; 47 34; 47 40; 47 42; 48 4; 48 8; 48 11; 48 41; 49 12; 49 20; 49 21; 49 36; 49 38; 49 42]
global d_x = [3.0, 5.0, 1.0, 2.0, 1.0, 3.0, 8.0, 9.0, 7.0, 6.0, 10.0, 1.0, 1.0, 6.0, 3.0, 1.0, 4.0, 10.0, 7.0, 6.0, 3.0, 3.0, 8.0, 5.0, 6.0, 7.0, 3.0, 7.0, 7.0, 9.0, 6.0, 6.0, 6.0, 3.0, 9.0, 6.0, 3.0, 8.0, 9.0, 9.0, 2.0, 2.0, 6.0, 1.0, 1.0, 2.0, 10.0, 2.0, 1.0, 10.0, 7.0, 5.0, 5.0, 9.0, 9.0, 1.0, 3.0, 6.0, 10.0, 3.0, 6.0, 3.0, 1.0, 4.0, 7.0, 1.0, 10.0, 8.0, 9.0, 7.0, 6.0, 8.0, 6.0, 3.0, 8.0, 5.0, 9.0, 7.0, 4.0, 2.0, 7.0, 10.0, 6.0, 5.0, 4.0, 8.0, 10.0, 7.0, 7.0, 9.0, 7.0, 8.0, 9.0, 6.0, 6.0, 9.0, 10.0, 9.0, 4.0, 6.0, 7.0, 10.0, 6.0, 7.0, 3.0, 1.0, 5.0, 3.0, 6.0, 2.0, 7.0, 1.0, 2.0, 7.0, 3.0, 6.0, 8.0, 6.0, 9.0, 5.0, 5.0, 10.0, 9.0, 1.0, 2.0, 8.0, 10.0, 6.0, 2.0, 2.0, 6.0, 10.0, 9.0, 2.0, 7.0, 4.0, 9.0, 7.0, 8.0, 8.0, 3.0, 9.0, 9.0, 4.0, 2.0, 7.0, 8.0, 3.0, 5.0, 9.0, 5.0, 6.0, 3.0, 4.0, 10.0, 2.0, 6.0, 2.0, 4.0, 8.0, 2.0, 8.0, 3.0, 8.0, 2.0, 4.0, 1.0, 8.0, 7.0, 1.0, 4.0, 2.0, 7.0, 7.0, 10.0, 3.0, 6.0, 2.0, 6.0, 9.0, 7.0, 1.0, 8.0, 4.0, 4.0, 7.0, 8.0, 5.0, 1.0, 6.0, 1.0, 1.0, 3.0, 7.0, 6.0, 10.0, 8.0, 4.0, 4.0, 9.0, 9.0, 9.0, 5.0, 6.0, 6.0, 9.0, 9.0, 9.0, 1.0, 1.0, 5.0, 4.0, 1.0, 5.0, 3.0, 2.0, 2.0, 2.0, 4.0, 9.0, 4.0, 7.0, 1.0, 3.0, 8.0, 1.0, 5.0, 4.0, 9.0, 5.0, 3.0, 5.0, 9.0, 6.0, 5.0, 7.0, 9.0, 8.0, 6.0, 5.0, 7.0, 2.0, 3.0, 5.0, 2.0, 1.0, 6.0, 7.0, 3.0, 6.0, 1.0, 6.0, 10.0, 4.0, 9.0, 10.0, 8.0, 7.0, 9.0, 5.0, 9.0, 10.0, 5.0, 3.0]
global b_x = 5
global d_y = [8.0, 2.0, 10.0, 4.0, 1.0, 5.0, 1.0, 6.0, 5.0, 3.0, 3.0, 7.0, 6.0, 3.0, 4.0, 5.0, 9.0, 1.0, 2.0, 7.0, 10.0, 2.0, 1.0, 6.0, 3.0, 2.0, 3.0, 7.0, 8.0, 6.0, 1.0, 3.0, 5.0, 5.0, 7.0, 9.0, 1.0, 5.0, 2.0, 7.0, 4.0, 2.0, 5.0, 2.0, 1.0, 3.0, 8.0, 6.0, 1.0, 5.0, 7.0, 8.0, 1.0, 4.0, 9.0, 9.0, 6.0, 2.0, 7.0, 1.0, 8.0, 7.0, 4.0, 10.0, 3.0, 10.0, 4.0, 2.0, 2.0, 5.0, 5.0, 6.0, 7.0, 8.0, 2.0, 3.0, 7.0, 6.0, 7.0, 8.0, 9.0, 4.0, 4.0, 8.0, 7.0, 4.0, 6.0, 6.0, 2.0, 8.0, 2.0, 2.0, 4.0, 10.0, 10.0, 7.0, 3.0, 10.0, 2.0, 5.0, 8.0, 8.0, 2.0, 2.0, 1.0, 1.0, 6.0, 4.0, 5.0, 3.0, 5.0, 10.0, 2.0, 3.0, 6.0, 7.0, 3.0, 5.0, 8.0, 10.0, 4.0, 6.0, 10.0, 9.0, 2.0, 4.0, 2.0, 3.0, 10.0, 6.0, 9.0, 5.0, 6.0, 5.0, 4.0, 3.0, 5.0, 6.0, 5.0, 1.0, 8.0, 1.0, 8.0, 7.0, 4.0, 5.0, 7.0, 7.0, 7.0, 2.0, 7.0, 8.0, 10.0, 5.0, 6.0, 5.0, 2.0, 6.0, 10.0, 6.0, 10.0, 2.0, 1.0, 7.0, 5.0, 7.0, 10.0, 3.0, 1.0, 10.0, 4.0, 5.0, 5.0, 6.0, 3.0, 4.0, 8.0, 2.0, 5.0, 10.0, 5.0, 7.0, 3.0, 8.0, 10.0, 7.0, 4.0, 10.0, 7.0, 2.0, 4.0, 7.0, 8.0, 7.0, 10.0, 5.0, 5.0, 5.0, 3.0, 7.0, 8.0, 4.0, 2.0, 7.0, 9.0, 6.0, 8.0, 4.0, 8.0, 9.0, 1.0, 4.0, 5.0, 3.0, 8.0, 1.0, 1.0, 2.0, 9.0, 6.0, 10.0, 2.0, 5.0, 4.0, 4.0, 7.0, 2.0, 4.0, 2.0, 8.0, 9.0, 4.0, 5.0, 1.0, 4.0, 6.0, 9.0, 9.0, 4.0, 10.0, 3.0, 3.0, 5.0, 2.0, 6.0, 3.0, 3.0, 5.0, 5.0, 5.0, 6.0, 8.0, 7.0, 7.0, 1.0, 9.0, 1.0, 3.0, 6.0, 5.0, 1.0, 8.0, 5.0, 10.0]
global b_y = 10
global p = [0.855, 0.032, 0.852, 0.186, 0.118, 0.769, 0.637, 0.469, 0.076, 0.827, 0.991, 0.061, 0.138, 0.897, 0.076, 0.704, 0.848, 0.248, 0.7, 0.914, 0.764, 0.818, 0.596, 0.068, 0.254, 0.969, 0.016, 0.817, 0.459, 0.35, 0.735, 0.747, 0.602, 0.679, 0.472, 0.393, 0.55, 0.064, 0.417, 0.911, 0.962, 0.464, 0.058, 0.981, 0.103, 0.202, 0.193, 0.49, 0.291, 0.355, 0.115, 0.686, 0.375, 0.325, 0.857, 0.9, 0.717, 0.811, 0.191, 0.889, 0.573, 0.063, 0.049, 0.288, 0.157, 0.265, 0.502, 0.308, 0.987, 0.425, 0.803, 0.864, 0.806, 0.562, 0.428, 0.799, 0.749, 0.232, 0.005, 0.533, 0.804, 0.888, 0.755, 0.057, 0.631, 0.654, 0.591, 0.4, 0.404, 0.433, 0.626, 0.001, 0.468, 0.099, 0.416, 0.289, 0.198, 0.714, 0.643, 0.111, 0.251, 0.283, 0.909, 0.497, 0.431, 0.774, 0.063, 0.387, 0.671, 0.844, 0.729, 0.286, 0.911, 0.047, 0.163, 0.569, 0.851, 0.788, 0.469, 0.929, 0.132, 0.514, 0.421, 0.81, 0.335, 0.29, 0.926, 0.901, 0.688, 0.954, 0.417, 0.528, 0.893, 0.673, 0.554, 0.45, 0.641, 0.309, 0.415, 0.921, 0.588, 0.751, 0.793, 0.651, 0.843, 0.669, 0.083, 0.64, 0.616, 0.736, 0.197, 0.187, 0.616, 0.221, 0.832, 0.68, 0.856, 0.119, 0.393, 0.677, 0.737, 0.366, 0.427, 0.47, 0.692, 0.407, 0.933, 0.812, 0.086, 0.364, 0.151, 0.655, 0.241, 0.198, 0.03, 0.554, 0.179, 0.45, 0.058, 0.858, 0.1, 0.802, 0.186, 0.243, 0.807, 0.951, 0.862, 0.946, 0.921, 0.78, 0.668, 0.311, 0.837, 0.208, 0.701, 0.677, 0.312, 0.881, 0.339, 0.934, 0.138, 0.154, 0.834, 0.485, 0.668, 0.83, 0.917, 0.344, 0.764, 0.862, 0.775, 0.613, 0.322, 0.854, 0.863, 0.164, 0.841, 0.427, 0.83, 0.587, 0.608, 0.7, 0.407, 0.94, 0.217, 0.818, 0.549, 0.275, 0.053, 0.548, 0.351, 0.24, 0.914, 0.887, 0.214, 0.377, 0.925, 0.153, 0.083, 0.826, 0.545, 0.985, 0.464, 0.52, 0.914, 0.501, 0.337, 0.418, 0.045, 0.386, 0.362, 0.473, 0.048, 0.721, 0.624, 0.191, 0.796, 0.557, 0.218, 0.862, 0.773, 0.85, 0.349, 0.995]
global q = [0.959, 0.863, 0.962, 0.488, 0.854, 0.912, 0.795, 0.51, 0.586, 0.847, 0.992, 0.963, 0.525, 0.914, 0.911, 0.936, 0.864, 0.954, 0.877, 0.934, 0.802, 0.849, 0.787, 0.669, 0.568, 0.986, 0.367, 0.936, 0.974, 0.4, 0.946, 0.983, 0.756, 0.991, 0.767, 0.471, 0.563, 0.111, 0.84, 0.996, 0.967, 0.769, 0.259, 0.999, 0.259, 0.807, 0.488, 0.835, 0.516, 0.94, 0.443, 0.871, 0.414, 0.704, 0.957, 0.904, 0.843, 0.938, 0.62, 0.985, 0.94, 0.342, 0.086, 0.809, 0.708, 0.32, 0.648, 0.665, 0.987, 0.946, 0.854, 0.939, 0.931, 0.934, 0.631, 0.977, 0.912, 0.419, 0.246, 0.939, 0.965, 0.944, 0.985, 0.77, 0.83, 0.932, 0.632, 0.954, 0.719, 0.502, 0.672, 0.316, 0.47, 0.587, 0.589, 0.746, 0.535, 0.816, 0.864, 0.156, 0.294, 0.362, 0.98, 0.778, 0.547, 0.942, 0.6, 0.693, 0.798, 0.893, 0.837, 0.729, 0.944, 0.233, 0.219, 0.896, 0.917, 0.975, 0.916, 0.965, 0.461, 0.648, 0.65, 0.878, 0.947, 0.647, 0.934, 0.942, 0.768, 0.969, 0.943, 0.917, 0.946, 0.705, 0.669, 0.998, 0.649, 0.694, 0.924, 0.968, 0.637, 0.938, 0.84, 0.721, 0.866, 0.972, 0.761, 0.758, 0.81, 0.916, 0.91, 0.23, 0.872, 0.445, 0.919, 0.947, 0.87, 0.646, 0.911, 0.956, 0.998, 0.497, 0.582, 0.567, 0.916, 0.768, 0.94, 0.92, 0.607, 0.448, 0.861, 0.683, 0.627, 0.692, 0.125, 0.709, 0.32, 0.454, 0.895, 0.91, 0.978, 0.874, 0.703, 0.637, 0.929, 0.955, 0.908, 0.996, 0.956, 0.785, 0.714, 0.551, 0.934, 0.59, 0.805, 0.702, 0.692, 0.938, 0.748, 0.983, 0.354, 0.547, 0.988, 0.903, 0.68, 0.883, 0.938, 0.72, 0.811, 0.901, 0.788, 0.683, 0.966, 0.896, 0.991, 0.571, 0.972, 0.713, 0.844, 0.865, 0.986, 0.826, 0.701, 0.978, 0.239, 0.834, 0.632, 0.981, 0.54, 0.851, 0.8, 0.337, 0.959, 0.984, 0.326, 0.439, 0.999, 0.733, 0.333, 0.967, 0.855, 0.991, 0.603, 0.82, 0.914, 0.701, 0.381, 0.699, 0.122, 0.731, 0.694, 0.888, 0.556, 0.903, 0.64, 0.391, 0.988, 0.584, 0.962, 0.973, 0.902, 0.992, 0.529, 0.996]
global origin = 1
global destination = 50