global arcs = [1 8; 1 13; 1 16; 1 18; 1 23; 1 49; 1 55; 1 57; 2 8; 2 18; 2 22; 2 39; 2 53; 2 57; 2 60; 3 12; 3 22; 3 23; 3 33; 3 57; 4 5; 4 22; 4 26; 4 32; 4 43; 4 44; 4 48; 4 52; 5 11; 5 18; 5 22; 5 26; 5 28; 5 32; 5 46; 5 47; 6 17; 6 20; 6 22; 6 35; 6 36; 6 39; 6 53; 6 55; 7 17; 7 30; 7 34; 7 37; 7 40; 7 53; 8 5; 8 10; 8 48; 9 2; 9 29; 9 30; 9 40; 9 44; 9 53; 9 56; 9 57; 10 6; 10 20; 10 22; 10 24; 10 29; 10 35; 10 55; 11 23; 11 26; 11 29; 11 40; 11 45; 11 46; 11 48; 12 5; 12 23; 12 56; 13 7; 13 8; 13 20; 13 29; 13 53; 13 56; 14 4; 14 10; 14 13; 14 15; 14 27; 14 32; 14 33; 15 2; 15 3; 15 26; 15 28; 15 31; 16 6; 16 24; 16 41; 16 49; 16 51; 16 54; 17 2; 17 11; 17 16; 17 25; 17 27; 17 28; 17 44; 17 49; 18 31; 18 34; 18 39; 19 8; 19 10; 19 27; 19 37; 19 45; 19 53; 20 3; 20 9; 20 12; 20 23; 20 29; 20 55; 21 4; 21 7; 21 12; 21 22; 21 25; 21 30; 21 36; 22 2; 22 6; 22 21; 22 23; 22 24; 22 27; 22 31; 22 37; 22 42; 22 60; 23 3; 23 4; 23 15; 23 20; 23 26; 23 41; 23 48; 24 7; 24 19; 24 23; 24 25; 24 27; 24 40; 24 41; 24 44; 24 51; 24 54; 25 9; 25 10; 25 22; 25 24; 25 34; 25 36; 25 39; 25 52; 26 6; 26 21; 26 28; 26 43; 27 22; 27 31; 27 53; 27 58; 28 4; 28 24; 28 27; 28 51; 29 9; 29 17; 29 48; 29 50; 30 23; 30 31; 30 52; 30 60; 31 14; 31 27; 31 35; 31 57; 32 8; 32 24; 32 34; 32 45; 32 48; 32 57; 33 9; 33 19; 33 35; 33 51; 33 52; 33 54; 34 23; 34 27; 34 29; 34 41; 34 49; 34 58; 35 23; 35 29; 35 32; 36 4; 36 22; 36 25; 36 54; 36 59; 37 6; 37 10; 37 20; 38 5; 38 6; 38 7; 38 24; 38 41; 38 42; 38 53; 39 16; 39 45; 39 49; 39 51; 40 16; 40 17; 41 4; 41 45; 41 52; 41 59; 42 3; 42 5; 42 29; 42 31; 43 8; 43 24; 43 37; 43 50; 44 24; 45 18; 45 22; 45 29; 45 32; 45 46; 45 56; 46 6; 46 13; 46 15; 46 18; 46 21; 46 23; 46 24; 46 39; 46 42; 46 49; 46 50; 46 54; 47 12; 47 29; 47 34; 47 42; 47 46; 47 50; 47 51; 48 7; 48 36; 49 2; 49 12; 49 16; 49 20; 49 34; 49 36; 49 37; 49 40; 49 53; 49 57; 50 9; 50 14; 50 36; 50 45; 50 48; 51 4; 51 25; 51 28; 51 34; 51 39; 52 15; 52 18; 52 19; 52 38; 52 39; 53 35; 53 55; 54 23; 54 55; 54 58; 55 2; 55 7; 55 11; 55 19; 55 22; 55 36; 55 40; 55 42; 55 54; 55 58; 56 7; 56 18; 56 21; 56 25; 56 44; 57 3; 57 38; 57 47; 57 52; 57 55; 57 60; 58 13; 58 42; 58 46; 58 47; 58 59; 59 4; 59 13; 59 15; 59 26; 59 27]
global d_x = [2.0, 4.0, 10.0, 10.0, 8.0, 8.0, 8.0, 2.0, 4.0, 9.0, 3.0, 3.0, 7.0, 6.0, 10.0, 5.0, 3.0, 7.0, 10.0, 10.0, 3.0, 4.0, 3.0, 2.0, 2.0, 2.0, 6.0, 3.0, 7.0, 6.0, 4.0, 9.0, 1.0, 7.0, 6.0, 5.0, 2.0, 6.0, 5.0, 9.0, 10.0, 4.0, 1.0, 2.0, 4.0, 4.0, 6.0, 10.0, 4.0, 6.0, 8.0, 5.0, 4.0, 10.0, 4.0, 4.0, 5.0, 5.0, 9.0, 2.0, 4.0, 5.0, 7.0, 9.0, 3.0, 10.0, 5.0, 8.0, 7.0, 8.0, 3.0, 8.0, 2.0, 6.0, 5.0, 5.0, 2.0, 5.0, 6.0, 5.0, 2.0, 7.0, 9.0, 5.0, 10.0, 3.0, 9.0, 6.0, 8.0, 1.0, 5.0, 5.0, 1.0, 10.0, 10.0, 9.0, 4.0, 1.0, 9.0, 10.0, 10.0, 7.0, 2.0, 4.0, 1.0, 9.0, 1.0, 10.0, 8.0, 10.0, 7.0, 3.0, 4.0, 8.0, 7.0, 5.0, 5.0, 3.0, 9.0, 6.0, 2.0, 5.0, 2.0, 4.0, 10.0, 3.0, 3.0, 9.0, 2.0, 8.0, 6.0, 9.0, 8.0, 3.0, 7.0, 8.0, 7.0, 8.0, 5.0, 8.0, 2.0, 6.0, 6.0, 9.0, 2.0, 7.0, 9.0, 7.0, 3.0, 1.0, 1.0, 4.0, 2.0, 10.0, 1.0, 2.0, 3.0, 9.0, 5.0, 2.0, 6.0, 3.0, 9.0, 7.0, 2.0, 8.0, 1.0, 7.0, 6.0, 4.0, 8.0, 3.0, 1.0, 1.0, 7.0, 9.0, 8.0, 4.0, 8.0, 3.0, 3.0, 6.0, 9.0, 2.0, 7.0, 1.0, 3.0, 5.0, 4.0, 9.0, 2.0, 7.0, 5.0, 6.0, 3.0, 9.0, 9.0, 5.0, 1.0, 5.0, 7.0, 3.0, 2.0, 5.0, 3.0, 6.0, 5.0, 6.0, 9.0, 9.0, 7.0, 6.0, 1.0, 9.0, 3.0, 7.0, 9.0, 10.0, 2.0, 1.0, 10.0, 7.0, 3.0, 1.0, 6.0, 4.0, 3.0, 5.0, 3.0, 1.0, 3.0, 6.0, 10.0, 1.0, 2.0, 2.0, 7.0, 9.0, 10.0, 3.0, 4.0, 2.0, 6.0, 10.0, 10.0, 1.0, 9.0, 6.0, 8.0, 6.0, 10.0, 4.0, 1.0, 10.0, 8.0, 4.0, 8.0, 9.0, 3.0, 4.0, 2.0, 8.0, 3.0, 10.0, 9.0, 4.0, 7.0, 5.0, 6.0, 5.0, 2.0, 3.0, 3.0, 4.0, 6.0, 1.0, 10.0, 9.0, 3.0, 1.0, 1.0, 2.0, 2.0, 3.0, 7.0, 10.0, 1.0, 5.0, 1.0, 5.0, 10.0, 5.0, 5.0, 4.0, 3.0, 10.0, 7.0, 4.0, 3.0, 10.0, 1.0, 3.0, 8.0, 1.0, 1.0, 9.0, 4.0, 5.0, 8.0, 7.0, 10.0, 5.0, 10.0, 8.0, 9.0, 7.0, 1.0, 1.0, 5.0, 9.0, 9.0, 5.0, 10.0, 2.0, 9.0, 1.0, 2.0, 6.0, 6.0, 5.0, 7.0, 5.0, 10.0, 1.0]
global b_x = 5
global d_y = [5.0, 2.0, 8.0, 1.0, 6.0, 2.0, 5.0, 2.0, 1.0, 7.0, 6.0, 10.0, 3.0, 4.0, 10.0, 4.0, 7.0, 9.0, 5.0, 7.0, 8.0, 8.0, 5.0, 8.0, 7.0, 10.0, 7.0, 4.0, 4.0, 3.0, 6.0, 1.0, 2.0, 9.0, 3.0, 4.0, 7.0, 10.0, 1.0, 4.0, 3.0, 7.0, 6.0, 5.0, 4.0, 2.0, 8.0, 10.0, 6.0, 10.0, 9.0, 7.0, 7.0, 1.0, 2.0, 4.0, 4.0, 8.0, 5.0, 1.0, 5.0, 1.0, 9.0, 1.0, 7.0, 4.0, 9.0, 3.0, 7.0, 10.0, 7.0, 3.0, 4.0, 8.0, 7.0, 3.0, 4.0, 8.0, 8.0, 3.0, 3.0, 5.0, 7.0, 8.0, 3.0, 10.0, 6.0, 3.0, 4.0, 10.0, 6.0, 7.0, 2.0, 8.0, 3.0, 6.0, 9.0, 1.0, 2.0, 10.0, 2.0, 2.0, 1.0, 9.0, 10.0, 2.0, 7.0, 1.0, 6.0, 8.0, 7.0, 1.0, 5.0, 9.0, 9.0, 10.0, 3.0, 9.0, 4.0, 7.0, 3.0, 7.0, 9.0, 6.0, 6.0, 2.0, 5.0, 8.0, 6.0, 8.0, 8.0, 6.0, 3.0, 9.0, 7.0, 5.0, 10.0, 2.0, 8.0, 3.0, 9.0, 8.0, 7.0, 5.0, 3.0, 3.0, 10.0, 10.0, 4.0, 3.0, 7.0, 9.0, 1.0, 5.0, 5.0, 6.0, 7.0, 6.0, 9.0, 9.0, 1.0, 10.0, 9.0, 5.0, 1.0, 10.0, 8.0, 2.0, 6.0, 7.0, 3.0, 6.0, 3.0, 8.0, 6.0, 8.0, 4.0, 8.0, 1.0, 10.0, 10.0, 5.0, 5.0, 4.0, 7.0, 7.0, 1.0, 10.0, 3.0, 1.0, 8.0, 5.0, 9.0, 6.0, 7.0, 5.0, 9.0, 10.0, 5.0, 2.0, 6.0, 2.0, 1.0, 4.0, 4.0, 4.0, 2.0, 6.0, 6.0, 9.0, 3.0, 7.0, 2.0, 10.0, 6.0, 10.0, 9.0, 2.0, 5.0, 3.0, 6.0, 5.0, 5.0, 3.0, 3.0, 4.0, 3.0, 2.0, 1.0, 5.0, 10.0, 5.0, 5.0, 8.0, 8.0, 6.0, 6.0, 8.0, 8.0, 1.0, 10.0, 5.0, 6.0, 6.0, 3.0, 2.0, 8.0, 8.0, 8.0, 3.0, 3.0, 4.0, 8.0, 3.0, 4.0, 3.0, 9.0, 6.0, 1.0, 1.0, 7.0, 7.0, 5.0, 9.0, 6.0, 1.0, 2.0, 3.0, 5.0, 1.0, 5.0, 3.0, 10.0, 2.0, 3.0, 5.0, 1.0, 1.0, 5.0, 3.0, 4.0, 3.0, 6.0, 7.0, 6.0, 1.0, 8.0, 7.0, 6.0, 7.0, 7.0, 10.0, 10.0, 1.0, 8.0, 1.0, 6.0, 7.0, 3.0, 2.0, 2.0, 5.0, 2.0, 2.0, 8.0, 1.0, 6.0, 8.0, 6.0, 6.0, 7.0, 2.0, 10.0, 10.0, 4.0, 9.0, 9.0, 3.0, 5.0, 3.0, 7.0, 1.0, 3.0, 4.0, 5.0, 1.0, 6.0, 3.0, 10.0, 7.0, 4.0, 4.0, 3.0, 8.0]
global b_y = 10
global p = [0.502, 0.655, 0.277, 0.815, 0.34, 0.591, 0.837, 0.005, 0.954, 0.583, 0.14, 0.015, 0.8, 0.76, 0.354, 0.43, 0.736, 0.392, 0.078, 0.379, 0.252, 0.271, 0.219, 0.153, 0.857, 0.41, 0.679, 0.507, 0.218, 0.018, 0.123, 0.613, 0.136, 0.576, 0.142, 0.097, 0.786, 0.773, 0.792, 0.557, 0.557, 0.889, 0.103, 0.033, 0.216, 0.362, 0.12, 0.371, 0.638, 0.26, 0.388, 0.577, 0.085, 0.43, 0.264, 0.06, 0.714, 0.97, 0.632, 0.862, 0.329, 0.253, 0.831, 0.954, 0.91, 0.545, 0.91, 0.872, 0.109, 0.005, 0.606, 0.007, 0.705, 0.359, 0.676, 0.251, 0.147, 0.105, 0.528, 0.938, 0.416, 0.229, 0.494, 0.591, 0.824, 0.547, 0.195, 0.481, 0.054, 0.036, 0.11, 0.634, 0.035, 0.341, 0.25, 0.284, 0.847, 0.371, 0.228, 0.505, 0.533, 0.803, 0.413, 0.725, 0.132, 0.825, 0.606, 0.242, 0.427, 0.47, 0.374, 0.952, 0.348, 0.089, 0.594, 0.096, 0.195, 0.585, 0.803, 0.183, 0.095, 0.257, 0.747, 0.874, 0.715, 0.977, 0.405, 0.029, 0.901, 0.053, 0.069, 0.657, 0.566, 0.871, 0.981, 0.221, 0.012, 0.858, 0.451, 0.104, 0.07, 0.833, 0.58, 0.329, 0.985, 0.543, 0.96, 0.601, 0.4, 0.91, 0.176, 0.599, 0.322, 0.755, 0.494, 0.213, 0.567, 0.399, 0.939, 0.719, 0.643, 0.012, 0.018, 0.572, 0.973, 0.506, 0.15, 0.454, 0.402, 0.106, 0.648, 0.376, 0.468, 0.063, 0.792, 0.556, 0.275, 0.869, 0.14, 0.927, 0.114, 0.349, 0.733, 0.413, 0.559, 0.261, 0.242, 0.492, 0.842, 0.021, 0.51, 0.37, 0.998, 0.605, 0.623, 0.188, 0.535, 0.008, 0.927, 0.697, 0.633, 0.49, 0.973, 0.732, 0.123, 0.035, 0.275, 0.267, 0.595, 0.895, 0.57, 0.278, 0.518, 0.721, 0.929, 0.665, 0.959, 0.131, 0.593, 0.947, 0.735, 0.873, 0.239, 0.677, 0.074, 0.57, 0.652, 0.711, 0.393, 0.891, 0.412, 0.577, 0.327, 0.933, 0.148, 0.031, 0.967, 0.377, 0.331, 0.746, 0.157, 0.061, 0.339, 0.835, 0.55, 0.917, 0.157, 0.03, 0.003, 0.429, 0.839, 0.644, 0.695, 0.552, 0.294, 0.627, 0.758, 0.415, 0.302, 0.019, 0.776, 0.651, 0.105, 0.894, 0.939, 0.023, 0.292, 0.841, 0.623, 0.092, 0.904, 0.543, 0.11, 0.183, 0.722, 0.187, 0.212, 0.289, 0.059, 0.701, 0.37, 0.712, 0.091, 0.716, 0.696, 0.34, 0.54, 0.26, 0.729, 0.579, 0.192, 0.904, 0.745, 0.662, 0.409, 0.485, 0.972, 0.017, 0.801, 0.937, 0.984, 0.944, 0.497, 0.876, 0.417, 0.888, 0.06, 0.81, 0.783, 0.294, 0.759, 0.914, 0.284, 0.563, 0.841, 0.018, 0.393, 0.947, 0.132, 0.482, 0.314, 0.102, 0.911, 0.187, 0.763, 0.393, 0.776, 0.104, 0.73, 0.995, 0.981, 0.192, 0.74, 0.299]
global q = [0.752, 0.742, 0.911, 0.989, 0.393, 0.978, 0.936, 0.791, 0.962, 0.602, 0.289, 0.677, 0.965, 0.914, 0.476, 0.496, 0.925, 0.685, 0.275, 0.744, 0.494, 0.316, 0.461, 0.862, 0.994, 0.797, 0.777, 0.547, 0.266, 0.018, 0.596, 0.715, 0.162, 0.955, 0.21, 0.685, 0.977, 0.949, 0.809, 0.921, 0.635, 0.933, 0.721, 0.744, 0.615, 0.869, 0.488, 0.556, 0.894, 0.656, 0.538, 0.963, 0.287, 0.932, 0.414, 0.197, 0.812, 0.99, 0.647, 0.862, 0.394, 0.567, 0.842, 0.997, 0.982, 0.63, 0.98, 0.909, 0.657, 0.668, 0.82, 0.304, 0.755, 0.442, 0.86, 0.37, 0.549, 0.273, 0.925, 0.99, 0.442, 0.259, 0.687, 0.747, 0.859, 0.676, 0.301, 0.499, 0.145, 0.457, 0.15, 0.805, 0.702, 0.342, 0.781, 0.419, 0.868, 0.816, 0.763, 0.706, 0.593, 0.914, 0.651, 0.876, 0.559, 0.958, 0.963, 0.287, 0.797, 0.525, 0.747, 0.988, 0.61, 0.689, 0.658, 0.915, 0.263, 0.806, 0.96, 0.769, 0.162, 0.699, 0.842, 0.988, 0.857, 0.99, 0.572, 0.174, 0.91, 0.085, 0.923, 0.909, 0.881, 0.921, 0.999, 0.231, 0.82, 0.945, 0.536, 0.194, 0.73, 0.937, 0.759, 0.591, 0.995, 0.815, 0.999, 0.77, 0.434, 0.939, 0.525, 0.657, 0.417, 0.8, 0.528, 0.614, 0.584, 0.615, 0.985, 0.773, 0.892, 0.754, 0.908, 0.778, 0.999, 0.799, 0.935, 0.773, 0.781, 0.132, 0.685, 0.853, 0.95, 0.74, 0.842, 0.584, 0.611, 0.941, 0.561, 0.946, 0.472, 0.722, 0.795, 0.926, 0.917, 0.494, 0.654, 0.549, 0.852, 0.727, 0.655, 0.86, 0.999, 0.864, 0.672, 0.901, 0.721, 0.916, 0.985, 0.763, 0.677, 0.793, 0.983, 0.787, 0.968, 0.686, 0.554, 0.748, 0.852, 0.897, 0.59, 0.357, 0.901, 0.74, 0.962, 0.786, 0.982, 0.563, 0.855, 0.989, 0.839, 0.916, 0.936, 0.988, 0.853, 0.932, 0.804, 0.971, 0.592, 0.908, 0.512, 0.964, 0.58, 0.946, 0.968, 0.575, 0.976, 0.708, 0.521, 0.813, 0.217, 0.315, 0.974, 0.99, 0.937, 0.948, 0.672, 0.574, 0.301, 0.841, 0.855, 0.769, 0.75, 0.621, 0.789, 0.769, 0.762, 0.799, 0.848, 0.511, 0.892, 0.889, 0.875, 0.921, 0.956, 0.334, 0.508, 0.93, 0.784, 0.528, 0.99, 0.752, 0.375, 0.56, 0.861, 0.216, 0.312, 0.393, 0.912, 0.879, 0.473, 0.993, 0.145, 0.925, 0.865, 0.405, 0.853, 0.509, 0.959, 0.74, 0.466, 0.946, 0.81, 0.762, 0.415, 0.992, 0.987, 0.314, 0.935, 0.946, 0.993, 0.994, 0.999, 0.958, 0.715, 0.995, 0.998, 0.966, 0.884, 0.902, 0.794, 0.933, 0.376, 0.809, 0.92, 0.033, 0.562, 0.987, 0.878, 0.716, 0.358, 0.421, 0.959, 0.888, 0.91, 0.98, 0.99, 0.955, 0.981, 0.999, 0.996, 0.885, 0.865, 0.704]
global origin = 1
global destination = 60