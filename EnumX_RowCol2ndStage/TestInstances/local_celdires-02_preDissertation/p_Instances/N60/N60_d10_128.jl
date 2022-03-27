global arcs = [1 4; 1 5; 1 40; 1 42; 1 57; 1 60; 2 24; 2 41; 3 18; 3 38; 3 42; 3 55; 3 56; 4 2; 4 8; 4 23; 4 31; 4 47; 5 15; 5 24; 5 37; 5 39; 5 49; 5 52; 5 55; 6 4; 6 12; 6 15; 6 53; 7 21; 7 31; 7 36; 8 10; 8 18; 8 38; 8 56; 8 59; 9 5; 9 6; 9 20; 9 24; 9 48; 10 3; 10 17; 10 27; 10 29; 10 34; 10 35; 10 40; 11 25; 11 40; 11 50; 12 4; 12 25; 12 49; 12 55; 13 17; 13 29; 13 30; 13 44; 13 50; 14 5; 14 15; 14 30; 15 40; 15 48; 16 2; 16 18; 16 27; 16 48; 16 50; 16 52; 17 13; 17 36; 17 50; 17 52; 17 59; 18 12; 18 24; 18 27; 18 31; 18 33; 18 40; 18 53; 19 2; 19 23; 19 30; 19 35; 19 47; 19 49; 20 3; 20 6; 20 14; 20 28; 20 35; 20 39; 20 47; 20 50; 20 53; 20 57; 20 59; 21 25; 21 26; 21 31; 21 53; 22 8; 22 11; 22 51; 23 3; 23 4; 23 34; 23 56; 24 9; 24 17; 24 35; 24 37; 24 54; 25 5; 25 29; 25 30; 25 56; 25 57; 25 58; 26 2; 26 16; 26 34; 26 38; 27 11; 27 14; 27 28; 27 40; 27 42; 27 48; 27 60; 28 3; 28 4; 28 44; 28 57; 28 59; 29 2; 29 4; 29 6; 29 16; 29 27; 29 28; 29 53; 30 6; 30 22; 30 34; 31 12; 31 14; 31 24; 31 29; 31 39; 31 45; 31 60; 32 2; 32 10; 32 17; 32 18; 32 60; 33 2; 33 15; 33 43; 33 48; 33 56; 33 60; 34 2; 34 3; 34 15; 34 16; 34 39; 34 40; 35 11; 35 24; 35 27; 35 28; 35 29; 35 30; 35 42; 35 46; 35 55; 36 5; 36 9; 36 44; 36 53; 36 56; 37 13; 37 16; 37 27; 37 30; 37 39; 38 22; 38 35; 38 44; 38 58; 38 60; 39 37; 39 46; 40 5; 40 8; 40 18; 40 24; 40 30; 40 51; 40 57; 41 3; 41 4; 41 5; 41 9; 41 13; 41 29; 41 34; 41 43; 41 47; 41 50; 41 53; 41 54; 42 6; 42 10; 42 12; 42 17; 42 27; 42 34; 42 47; 43 10; 43 30; 43 32; 43 39; 43 44; 44 3; 44 10; 44 22; 44 32; 45 22; 45 24; 45 25; 45 43; 45 51; 45 54; 46 5; 46 6; 46 8; 46 31; 46 45; 46 51; 47 3; 47 5; 47 14; 47 17; 47 27; 47 42; 48 16; 48 20; 48 49; 48 53; 48 54; 49 7; 49 13; 49 15; 49 24; 49 51; 50 5; 50 25; 50 34; 50 37; 50 39; 50 57; 51 7; 51 27; 51 43; 51 54; 52 21; 52 26; 52 47; 52 60; 53 10; 53 16; 53 18; 53 19; 53 42; 54 5; 54 16; 54 23; 54 35; 54 38; 54 41; 54 49; 54 55; 55 6; 55 33; 55 36; 55 38; 55 44; 55 52; 55 58; 56 4; 56 11; 56 27; 56 46; 56 54; 56 58; 57 2; 57 8; 57 18; 57 26; 57 33; 57 47; 58 8; 58 20; 58 28; 58 40; 58 43; 59 15; 59 33; 59 52]
global d_x = [7.0, 7.0, 1.0, 7.0, 1.0, 5.0, 10.0, 4.0, 7.0, 4.0, 4.0, 10.0, 7.0, 10.0, 4.0, 3.0, 3.0, 4.0, 4.0, 9.0, 6.0, 9.0, 5.0, 5.0, 6.0, 3.0, 9.0, 7.0, 9.0, 3.0, 1.0, 2.0, 8.0, 9.0, 2.0, 7.0, 4.0, 6.0, 5.0, 2.0, 6.0, 5.0, 2.0, 5.0, 9.0, 1.0, 6.0, 4.0, 9.0, 5.0, 1.0, 6.0, 9.0, 4.0, 6.0, 5.0, 1.0, 3.0, 9.0, 6.0, 6.0, 9.0, 6.0, 3.0, 2.0, 10.0, 5.0, 8.0, 1.0, 6.0, 2.0, 8.0, 3.0, 3.0, 7.0, 6.0, 9.0, 8.0, 7.0, 3.0, 5.0, 3.0, 1.0, 9.0, 5.0, 8.0, 7.0, 9.0, 5.0, 1.0, 2.0, 5.0, 9.0, 5.0, 6.0, 1.0, 10.0, 10.0, 6.0, 3.0, 8.0, 8.0, 1.0, 8.0, 10.0, 3.0, 4.0, 5.0, 7.0, 4.0, 3.0, 3.0, 9.0, 7.0, 7.0, 5.0, 3.0, 2.0, 5.0, 7.0, 6.0, 4.0, 3.0, 4.0, 8.0, 3.0, 10.0, 9.0, 2.0, 5.0, 5.0, 5.0, 4.0, 10.0, 10.0, 3.0, 9.0, 2.0, 2.0, 7.0, 8.0, 1.0, 5.0, 5.0, 10.0, 10.0, 4.0, 8.0, 4.0, 8.0, 7.0, 8.0, 4.0, 6.0, 5.0, 10.0, 1.0, 9.0, 5.0, 9.0, 4.0, 3.0, 3.0, 7.0, 9.0, 4.0, 6.0, 2.0, 6.0, 1.0, 1.0, 9.0, 2.0, 5.0, 7.0, 4.0, 1.0, 2.0, 3.0, 1.0, 3.0, 5.0, 4.0, 3.0, 9.0, 5.0, 3.0, 7.0, 4.0, 5.0, 7.0, 4.0, 2.0, 3.0, 10.0, 8.0, 6.0, 9.0, 10.0, 1.0, 5.0, 10.0, 8.0, 6.0, 6.0, 9.0, 10.0, 9.0, 5.0, 7.0, 7.0, 6.0, 6.0, 3.0, 8.0, 8.0, 2.0, 5.0, 3.0, 9.0, 6.0, 1.0, 6.0, 3.0, 6.0, 5.0, 10.0, 8.0, 1.0, 3.0, 8.0, 3.0, 5.0, 8.0, 10.0, 3.0, 6.0, 3.0, 1.0, 6.0, 8.0, 7.0, 8.0, 10.0, 1.0, 3.0, 3.0, 3.0, 4.0, 9.0, 4.0, 1.0, 9.0, 1.0, 5.0, 6.0, 9.0, 1.0, 4.0, 7.0, 10.0, 3.0, 4.0, 10.0, 10.0, 8.0, 6.0, 10.0, 1.0, 1.0, 6.0, 3.0, 6.0, 10.0, 8.0, 9.0, 7.0, 9.0, 3.0, 7.0, 4.0, 6.0, 3.0, 7.0, 5.0, 7.0, 5.0, 2.0, 5.0, 2.0, 1.0, 2.0, 5.0, 8.0, 2.0, 5.0, 5.0, 8.0, 6.0, 7.0, 4.0, 8.0, 9.0, 4.0, 1.0, 10.0, 7.0, 6.0, 6.0, 10.0, 9.0, 1.0, 7.0, 3.0, 3.0, 7.0]
global b_x = 5
global d_y = [4.0, 1.0, 4.0, 1.0, 6.0, 6.0, 8.0, 10.0, 4.0, 5.0, 1.0, 4.0, 7.0, 7.0, 3.0, 7.0, 6.0, 6.0, 4.0, 10.0, 5.0, 4.0, 2.0, 7.0, 10.0, 4.0, 9.0, 8.0, 4.0, 10.0, 5.0, 4.0, 9.0, 5.0, 5.0, 8.0, 6.0, 4.0, 2.0, 6.0, 8.0, 1.0, 8.0, 2.0, 7.0, 3.0, 8.0, 10.0, 10.0, 9.0, 6.0, 7.0, 2.0, 6.0, 2.0, 5.0, 2.0, 3.0, 6.0, 9.0, 6.0, 1.0, 4.0, 8.0, 9.0, 4.0, 3.0, 7.0, 9.0, 6.0, 6.0, 6.0, 9.0, 2.0, 7.0, 7.0, 3.0, 9.0, 2.0, 8.0, 3.0, 8.0, 9.0, 5.0, 5.0, 7.0, 1.0, 7.0, 10.0, 7.0, 7.0, 5.0, 6.0, 7.0, 10.0, 4.0, 4.0, 2.0, 9.0, 5.0, 4.0, 7.0, 10.0, 3.0, 10.0, 8.0, 2.0, 1.0, 5.0, 1.0, 1.0, 6.0, 6.0, 2.0, 10.0, 6.0, 10.0, 4.0, 2.0, 4.0, 10.0, 2.0, 9.0, 6.0, 8.0, 10.0, 7.0, 2.0, 3.0, 10.0, 3.0, 6.0, 3.0, 3.0, 1.0, 7.0, 4.0, 1.0, 1.0, 8.0, 5.0, 8.0, 1.0, 6.0, 5.0, 7.0, 5.0, 6.0, 10.0, 6.0, 8.0, 4.0, 9.0, 7.0, 10.0, 9.0, 3.0, 9.0, 7.0, 5.0, 1.0, 9.0, 10.0, 2.0, 5.0, 8.0, 2.0, 9.0, 5.0, 3.0, 5.0, 4.0, 9.0, 7.0, 1.0, 9.0, 10.0, 1.0, 8.0, 3.0, 2.0, 10.0, 6.0, 4.0, 3.0, 4.0, 8.0, 10.0, 5.0, 10.0, 8.0, 3.0, 2.0, 10.0, 9.0, 5.0, 10.0, 10.0, 6.0, 3.0, 5.0, 9.0, 9.0, 7.0, 6.0, 1.0, 1.0, 10.0, 2.0, 5.0, 3.0, 9.0, 10.0, 6.0, 9.0, 3.0, 4.0, 9.0, 8.0, 1.0, 4.0, 8.0, 4.0, 9.0, 5.0, 4.0, 1.0, 9.0, 9.0, 10.0, 3.0, 8.0, 2.0, 3.0, 6.0, 10.0, 8.0, 5.0, 6.0, 8.0, 5.0, 5.0, 3.0, 10.0, 4.0, 10.0, 8.0, 10.0, 7.0, 10.0, 3.0, 2.0, 9.0, 2.0, 4.0, 5.0, 7.0, 5.0, 5.0, 9.0, 7.0, 8.0, 3.0, 9.0, 6.0, 4.0, 5.0, 3.0, 1.0, 6.0, 8.0, 9.0, 3.0, 5.0, 5.0, 10.0, 5.0, 8.0, 8.0, 9.0, 4.0, 5.0, 5.0, 5.0, 9.0, 4.0, 6.0, 9.0, 4.0, 2.0, 5.0, 6.0, 6.0, 5.0, 8.0, 8.0, 2.0, 6.0, 10.0, 9.0, 7.0, 2.0, 5.0, 4.0, 2.0, 5.0, 7.0, 5.0, 8.0, 6.0, 6.0, 4.0, 3.0, 1.0, 10.0, 1.0]
global b_y = 10
global p = [0.825, 0.652, 0.641, 0.789, 0.191, 0.692, 0.622, 0.655, 0.095, 0.511, 0.522, 0.509, 0.288, 0.011, 0.099, 0.671, 0.897, 0.799, 0.954, 0.093, 0.041, 0.512, 0.99, 0.093, 0.298, 0.346, 0.669, 0.23, 0.436, 0.739, 0.392, 0.749, 0.071, 0.792, 0.265, 0.392, 0.712, 0.02, 0.522, 0.331, 0.092, 0.194, 0.216, 0.813, 0.582, 0.865, 0.692, 0.62, 0.414, 0.886, 0.536, 0.793, 0.079, 0.504, 0.186, 0.723, 0.836, 0.755, 0.935, 0.776, 0.297, 0.6, 0.014, 0.369, 0.35, 0.967, 0.298, 0.005, 0.586, 0.212, 0.313, 0.602, 0.136, 0.748, 0.433, 0.498, 0.652, 0.497, 0.211, 0.128, 0.087, 0.185, 0.096, 0.498, 0.096, 0.243, 0.456, 0.229, 0.39, 0.663, 0.007, 0.035, 0.378, 0.309, 0.349, 0.079, 0.145, 0.084, 0.356, 0.012, 0.34, 0.173, 0.732, 0.608, 0.682, 0.23, 0.16, 0.336, 0.197, 0.547, 0.446, 0.282, 0.415, 0.461, 0.592, 0.186, 0.925, 0.219, 0.168, 0.154, 0.179, 0.87, 0.815, 0.015, 0.534, 0.357, 0.428, 0.166, 0.845, 0.562, 0.564, 0.15, 0.13, 0.609, 0.273, 0.658, 0.806, 0.387, 0.356, 0.449, 0.423, 0.787, 0.443, 0.866, 0.945, 0.256, 0.933, 0.243, 0.943, 0.14, 0.468, 0.105, 0.899, 0.271, 0.278, 0.911, 0.177, 0.473, 0.472, 0.513, 0.271, 0.497, 0.049, 0.114, 0.878, 0.267, 0.556, 0.987, 0.091, 0.556, 0.041, 0.375, 0.622, 0.355, 0.515, 0.367, 0.784, 0.533, 0.186, 0.396, 0.271, 0.996, 0.556, 0.191, 0.395, 0.792, 0.368, 0.657, 0.232, 0.256, 0.367, 0.251, 0.983, 0.438, 0.985, 0.701, 0.013, 0.848, 0.734, 0.225, 0.567, 0.948, 0.044, 0.778, 0.236, 0.051, 0.11, 0.96, 0.345, 0.185, 0.396, 0.932, 0.385, 0.335, 0.456, 0.713, 0.677, 0.409, 0.52, 0.223, 0.604, 0.154, 0.499, 0.981, 0.827, 0.351, 0.999, 0.191, 0.553, 0.473, 0.098, 0.213, 0.326, 0.98, 0.548, 0.796, 0.539, 0.859, 0.522, 0.339, 0.868, 0.539, 0.694, 0.43, 0.742, 0.987, 0.05, 0.774, 0.114, 0.368, 0.502, 0.206, 0.302, 0.549, 0.613, 0.584, 0.573, 0.671, 0.763, 0.689, 0.627, 0.416, 0.511, 0.575, 0.426, 0.532, 0.183, 0.702, 0.218, 0.345, 0.87, 0.253, 0.197, 0.622, 0.156, 0.641, 0.786, 0.908, 0.06, 0.298, 0.074, 0.464, 0.489, 0.376, 0.866, 0.835, 0.394, 0.179, 0.342, 0.057, 0.852, 0.779, 0.571, 0.989, 0.089, 0.837, 0.288, 0.975, 0.143, 0.967, 0.772, 0.213, 0.193, 0.505, 0.749, 0.488, 0.905, 0.009, 0.831, 0.482, 0.173, 0.43, 0.874, 0.824, 0.615, 0.291]
global q = [0.985, 0.715, 0.722, 0.908, 0.528, 0.93, 0.73, 0.936, 0.841, 0.599, 0.939, 0.881, 0.679, 0.708, 0.966, 0.956, 0.898, 0.842, 0.992, 0.547, 0.193, 0.901, 0.997, 0.453, 0.921, 0.954, 0.815, 0.671, 0.543, 0.875, 0.65, 0.985, 0.821, 0.961, 0.35, 0.569, 0.803, 0.273, 0.668, 0.549, 0.707, 0.58, 0.513, 0.97, 0.95, 0.906, 0.873, 0.921, 0.7, 0.956, 0.878, 0.822, 0.802, 0.513, 0.441, 0.993, 0.885, 0.758, 0.956, 0.899, 0.383, 0.804, 0.797, 0.686, 0.721, 0.967, 0.732, 0.128, 0.647, 0.317, 0.945, 0.877, 0.5, 0.792, 0.435, 0.638, 0.847, 0.525, 0.88, 0.512, 0.462, 0.68, 0.355, 0.572, 0.267, 0.35, 0.727, 0.597, 0.639, 0.84, 0.455, 0.354, 0.664, 0.755, 0.538, 0.478, 0.37, 0.313, 0.469, 0.565, 0.785, 0.348, 0.985, 0.744, 0.886, 0.758, 0.466, 0.59, 0.503, 0.744, 0.702, 0.755, 0.564, 0.639, 0.613, 0.429, 0.998, 0.408, 0.291, 0.99, 0.622, 0.924, 0.9, 0.851, 0.746, 0.615, 0.796, 0.482, 0.904, 0.888, 0.662, 0.425, 0.333, 0.624, 0.395, 0.727, 0.894, 0.7, 0.668, 0.63, 0.857, 0.98, 0.523, 0.884, 0.95, 0.973, 0.984, 0.67, 0.994, 0.292, 0.648, 0.259, 0.956, 0.429, 0.294, 0.929, 0.547, 0.921, 0.793, 0.758, 0.458, 0.905, 0.676, 0.572, 0.913, 0.279, 0.837, 0.99, 0.408, 0.919, 0.591, 0.933, 0.968, 0.959, 0.968, 0.662, 0.925, 0.692, 0.716, 0.7, 0.444, 0.998, 0.733, 0.339, 0.869, 0.894, 0.888, 0.973, 0.998, 0.881, 0.739, 0.379, 0.99, 0.747, 0.988, 0.701, 0.022, 0.975, 0.895, 0.233, 0.854, 0.963, 0.677, 0.837, 0.613, 0.104, 0.467, 0.992, 0.552, 0.201, 0.849, 0.981, 0.993, 0.844, 0.505, 0.732, 0.916, 0.74, 0.776, 0.869, 0.832, 0.426, 0.542, 0.987, 0.836, 0.637, 0.999, 0.702, 0.648, 0.549, 0.971, 0.391, 0.555, 0.983, 0.945, 0.879, 0.808, 0.945, 0.652, 0.671, 0.879, 0.638, 0.78, 0.508, 0.857, 0.99, 0.492, 0.967, 0.383, 0.872, 0.902, 0.7, 0.545, 0.701, 0.87, 0.628, 0.802, 0.823, 0.829, 0.947, 0.667, 0.959, 0.899, 0.769, 0.489, 0.763, 0.292, 0.764, 0.552, 0.775, 0.912, 0.402, 0.862, 0.932, 0.398, 0.956, 0.908, 0.994, 0.898, 0.463, 0.393, 0.982, 0.765, 0.938, 0.991, 0.919, 0.767, 0.541, 0.888, 0.265, 0.958, 0.892, 0.755, 0.997, 0.837, 0.981, 0.507, 0.999, 0.842, 0.998, 0.956, 0.754, 0.421, 0.675, 0.911, 0.989, 0.943, 0.215, 0.953, 0.799, 0.864, 0.457, 0.929, 0.827, 0.783, 0.709]
global origin = 1
global destination = 60