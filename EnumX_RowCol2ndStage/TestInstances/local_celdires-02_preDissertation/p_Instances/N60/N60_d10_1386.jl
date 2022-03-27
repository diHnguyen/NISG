global arcs = [1 4; 1 22; 1 30; 1 38; 1 41; 1 46; 1 58; 2 6; 2 7; 2 10; 2 60; 3 6; 3 7; 3 24; 3 30; 3 38; 3 48; 3 58; 4 5; 4 14; 4 27; 5 2; 5 6; 5 11; 5 16; 5 19; 5 22; 5 28; 5 36; 5 38; 5 54; 6 3; 6 4; 6 31; 6 33; 6 40; 7 5; 7 10; 7 12; 7 32; 7 42; 7 45; 7 47; 7 52; 8 11; 8 13; 8 47; 8 55; 9 17; 9 18; 9 33; 10 12; 10 17; 10 22; 10 23; 10 24; 10 29; 10 32; 10 46; 10 51; 10 52; 10 55; 11 8; 11 34; 11 37; 11 40; 11 42; 11 44; 12 2; 12 4; 12 39; 12 41; 12 50; 12 53; 13 4; 13 9; 13 11; 13 20; 13 34; 14 36; 14 46; 14 49; 14 53; 15 6; 15 10; 15 18; 15 24; 15 46; 15 58; 16 45; 16 53; 17 4; 17 7; 17 18; 17 22; 17 23; 17 38; 18 12; 18 33; 18 39; 18 45; 18 57; 19 11; 19 32; 19 33; 19 46; 20 17; 20 34; 20 43; 20 52; 20 54; 21 10; 21 14; 21 32; 21 41; 21 49; 21 53; 22 5; 22 17; 22 23; 22 27; 22 58; 23 8; 23 29; 23 33; 23 35; 23 36; 23 38; 23 49; 24 4; 24 8; 24 10; 24 22; 24 48; 24 57; 25 24; 25 26; 25 29; 25 51; 25 53; 25 58; 26 12; 26 31; 26 41; 27 4; 27 17; 27 21; 27 30; 27 38; 28 2; 28 3; 28 9; 28 52; 28 53; 29 14; 29 19; 29 27; 29 32; 29 44; 30 2; 30 8; 30 27; 30 35; 30 44; 30 45; 31 6; 31 16; 31 20; 32 10; 32 12; 32 30; 32 40; 32 49; 33 15; 33 19; 33 34; 33 39; 33 46; 33 48; 33 52; 34 7; 34 15; 34 18; 34 21; 34 23; 34 25; 34 37; 34 38; 35 5; 35 13; 35 22; 35 38; 35 48; 35 51; 36 9; 36 23; 36 33; 36 38; 37 3; 37 13; 37 18; 37 29; 37 32; 37 44; 37 53; 38 2; 38 12; 38 14; 38 32; 38 58; 38 60; 39 8; 39 54; 40 9; 40 18; 40 31; 40 35; 40 37; 40 41; 40 48; 40 53; 41 4; 41 6; 41 18; 41 22; 41 33; 41 44; 41 51; 41 53; 42 16; 42 22; 42 26; 42 34; 42 46; 42 50; 42 51; 42 52; 42 53; 43 4; 43 8; 43 9; 43 12; 43 27; 43 30; 43 32; 43 37; 43 52; 44 33; 44 43; 44 47; 45 3; 45 10; 45 15; 45 20; 45 21; 45 30; 45 47; 46 7; 46 29; 46 36; 46 41; 47 11; 47 43; 47 55; 48 14; 48 43; 48 57; 49 11; 49 24; 49 56; 50 37; 50 41; 50 57; 51 3; 51 6; 51 9; 51 30; 51 36; 51 40; 52 4; 52 22; 52 27; 53 4; 53 11; 53 19; 53 21; 53 31; 54 3; 54 5; 54 17; 54 22; 54 34; 54 41; 54 59; 55 18; 55 19; 55 25; 55 32; 55 54; 56 8; 56 9; 56 17; 56 40; 56 47; 56 55; 57 2; 57 18; 57 22; 57 44; 57 45; 57 49; 57 54; 58 5; 58 13; 58 27; 58 34; 58 35; 58 42; 58 43; 58 53; 58 55; 58 60; 59 13; 59 21; 59 25; 59 33; 59 51; 59 57]
global d_x = [4.0, 7.0, 2.0, 5.0, 6.0, 6.0, 8.0, 3.0, 1.0, 8.0, 9.0, 5.0, 4.0, 6.0, 3.0, 7.0, 7.0, 1.0, 10.0, 6.0, 9.0, 8.0, 4.0, 9.0, 5.0, 3.0, 8.0, 2.0, 8.0, 1.0, 8.0, 9.0, 5.0, 4.0, 3.0, 1.0, 9.0, 5.0, 3.0, 4.0, 1.0, 1.0, 7.0, 1.0, 8.0, 4.0, 6.0, 2.0, 1.0, 4.0, 2.0, 8.0, 4.0, 6.0, 2.0, 6.0, 6.0, 6.0, 10.0, 10.0, 6.0, 3.0, 9.0, 6.0, 3.0, 9.0, 7.0, 8.0, 6.0, 8.0, 9.0, 5.0, 9.0, 1.0, 7.0, 5.0, 5.0, 5.0, 7.0, 5.0, 9.0, 2.0, 3.0, 2.0, 7.0, 2.0, 7.0, 10.0, 5.0, 5.0, 9.0, 8.0, 7.0, 3.0, 10.0, 10.0, 10.0, 2.0, 7.0, 8.0, 3.0, 3.0, 2.0, 8.0, 1.0, 7.0, 7.0, 3.0, 7.0, 3.0, 4.0, 1.0, 4.0, 7.0, 3.0, 9.0, 4.0, 1.0, 2.0, 9.0, 10.0, 9.0, 9.0, 8.0, 6.0, 8.0, 3.0, 4.0, 6.0, 1.0, 5.0, 1.0, 5.0, 6.0, 5.0, 6.0, 2.0, 2.0, 5.0, 5.0, 5.0, 7.0, 9.0, 6.0, 7.0, 8.0, 10.0, 1.0, 6.0, 1.0, 2.0, 3.0, 5.0, 4.0, 7.0, 9.0, 10.0, 5.0, 10.0, 7.0, 1.0, 10.0, 7.0, 10.0, 1.0, 6.0, 2.0, 8.0, 6.0, 3.0, 1.0, 5.0, 5.0, 10.0, 8.0, 7.0, 1.0, 6.0, 6.0, 4.0, 5.0, 5.0, 5.0, 10.0, 9.0, 2.0, 6.0, 2.0, 1.0, 5.0, 9.0, 10.0, 5.0, 4.0, 2.0, 10.0, 7.0, 4.0, 7.0, 4.0, 9.0, 4.0, 6.0, 8.0, 8.0, 1.0, 5.0, 9.0, 1.0, 8.0, 8.0, 1.0, 5.0, 10.0, 7.0, 2.0, 6.0, 6.0, 3.0, 1.0, 4.0, 7.0, 7.0, 1.0, 1.0, 8.0, 2.0, 9.0, 2.0, 10.0, 10.0, 3.0, 1.0, 8.0, 5.0, 6.0, 9.0, 2.0, 6.0, 10.0, 4.0, 7.0, 7.0, 3.0, 5.0, 4.0, 3.0, 5.0, 10.0, 1.0, 1.0, 6.0, 7.0, 3.0, 10.0, 5.0, 4.0, 2.0, 4.0, 3.0, 1.0, 5.0, 8.0, 5.0, 5.0, 6.0, 3.0, 4.0, 7.0, 7.0, 2.0, 2.0, 10.0, 5.0, 5.0, 2.0, 4.0, 5.0, 1.0, 3.0, 10.0, 10.0, 4.0, 6.0, 2.0, 4.0, 8.0, 5.0, 5.0, 4.0, 4.0, 5.0, 8.0, 3.0, 5.0, 3.0, 5.0, 8.0, 2.0, 1.0, 1.0, 10.0, 9.0, 5.0, 10.0, 10.0, 8.0, 3.0, 8.0, 1.0, 1.0, 6.0, 7.0, 8.0, 9.0, 8.0, 8.0, 7.0, 8.0, 9.0, 2.0, 3.0, 8.0, 4.0, 9.0, 3.0, 4.0, 2.0]
global b_x = 5
global d_y = [5.0, 9.0, 8.0, 6.0, 4.0, 8.0, 4.0, 4.0, 1.0, 5.0, 4.0, 2.0, 3.0, 6.0, 1.0, 2.0, 7.0, 7.0, 10.0, 5.0, 10.0, 4.0, 10.0, 10.0, 1.0, 10.0, 1.0, 1.0, 1.0, 4.0, 3.0, 2.0, 6.0, 5.0, 4.0, 2.0, 10.0, 9.0, 2.0, 9.0, 7.0, 5.0, 1.0, 6.0, 3.0, 3.0, 7.0, 9.0, 8.0, 10.0, 4.0, 6.0, 7.0, 6.0, 9.0, 6.0, 4.0, 4.0, 8.0, 6.0, 8.0, 1.0, 9.0, 4.0, 7.0, 6.0, 5.0, 4.0, 10.0, 2.0, 4.0, 7.0, 7.0, 10.0, 1.0, 9.0, 4.0, 4.0, 7.0, 2.0, 9.0, 1.0, 8.0, 5.0, 9.0, 8.0, 10.0, 10.0, 2.0, 8.0, 6.0, 7.0, 3.0, 6.0, 10.0, 3.0, 3.0, 3.0, 5.0, 6.0, 10.0, 4.0, 1.0, 3.0, 5.0, 10.0, 4.0, 3.0, 2.0, 9.0, 6.0, 8.0, 9.0, 7.0, 2.0, 5.0, 9.0, 8.0, 4.0, 4.0, 9.0, 4.0, 1.0, 2.0, 6.0, 7.0, 2.0, 8.0, 7.0, 7.0, 9.0, 5.0, 2.0, 6.0, 4.0, 1.0, 5.0, 2.0, 1.0, 6.0, 8.0, 2.0, 1.0, 8.0, 6.0, 2.0, 6.0, 2.0, 10.0, 3.0, 3.0, 8.0, 10.0, 8.0, 6.0, 10.0, 9.0, 2.0, 1.0, 2.0, 9.0, 5.0, 4.0, 8.0, 2.0, 5.0, 5.0, 5.0, 4.0, 2.0, 9.0, 2.0, 8.0, 5.0, 10.0, 9.0, 9.0, 3.0, 8.0, 4.0, 10.0, 6.0, 10.0, 7.0, 7.0, 3.0, 3.0, 8.0, 2.0, 4.0, 2.0, 9.0, 9.0, 7.0, 3.0, 7.0, 9.0, 4.0, 1.0, 1.0, 4.0, 7.0, 8.0, 7.0, 6.0, 3.0, 9.0, 8.0, 5.0, 9.0, 8.0, 6.0, 9.0, 1.0, 5.0, 6.0, 10.0, 4.0, 10.0, 3.0, 1.0, 4.0, 5.0, 3.0, 3.0, 8.0, 4.0, 10.0, 7.0, 10.0, 9.0, 7.0, 8.0, 3.0, 3.0, 3.0, 8.0, 2.0, 8.0, 8.0, 4.0, 8.0, 5.0, 1.0, 3.0, 3.0, 3.0, 10.0, 3.0, 1.0, 2.0, 7.0, 5.0, 3.0, 5.0, 10.0, 2.0, 9.0, 1.0, 7.0, 8.0, 6.0, 2.0, 10.0, 1.0, 4.0, 5.0, 9.0, 6.0, 3.0, 6.0, 2.0, 4.0, 2.0, 8.0, 8.0, 4.0, 6.0, 8.0, 7.0, 6.0, 2.0, 9.0, 5.0, 5.0, 5.0, 3.0, 5.0, 4.0, 5.0, 9.0, 2.0, 8.0, 2.0, 2.0, 1.0, 3.0, 9.0, 1.0, 3.0, 2.0, 5.0, 8.0, 2.0, 1.0, 7.0, 5.0, 7.0, 8.0, 10.0, 6.0, 6.0, 1.0, 7.0, 3.0, 4.0, 3.0, 4.0, 4.0, 5.0, 8.0, 9.0, 5.0, 7.0, 4.0, 1.0, 3.0, 6.0]
global b_y = 10
global p = [0.007, 0.23, 0.506, 0.027, 0.883, 0.654, 0.023, 0.724, 0.266, 0.822, 0.805, 0.555, 0.553, 0.635, 0.696, 0.421, 0.439, 0.438, 0.549, 0.012, 0.275, 0.97, 0.402, 0.874, 0.783, 0.261, 0.402, 0.207, 0.284, 0.278, 0.135, 0.216, 0.32, 0.874, 0.824, 0.058, 0.001, 0.352, 0.471, 0.18, 0.254, 0.98, 0.651, 0.626, 0.12, 0.922, 0.463, 0.083, 0.036, 0.995, 0.195, 0.932, 0.629, 0.913, 0.161, 0.485, 0.868, 0.957, 0.483, 0.181, 0.245, 0.621, 0.697, 0.473, 0.631, 0.66, 0.652, 0.139, 0.753, 0.869, 0.076, 0.153, 0.329, 0.784, 0.678, 0.015, 0.662, 0.582, 0.991, 0.704, 0.283, 0.274, 0.289, 0.524, 0.708, 0.915, 0.065, 0.926, 0.11, 0.744, 0.215, 0.025, 0.901, 0.614, 0.236, 0.994, 0.679, 0.133, 0.247, 0.34, 0.118, 0.649, 0.994, 0.42, 0.956, 0.542, 0.777, 0.59, 0.935, 0.198, 0.955, 0.154, 0.468, 0.774, 0.955, 0.456, 0.156, 0.101, 0.146, 0.155, 0.462, 0.278, 0.348, 0.602, 0.496, 0.645, 0.529, 0.238, 0.885, 0.375, 0.251, 0.806, 0.963, 0.099, 0.724, 0.634, 0.64, 0.396, 0.456, 0.445, 0.093, 0.843, 0.654, 0.846, 0.072, 0.315, 0.516, 0.029, 0.118, 0.67, 0.957, 0.776, 0.228, 0.252, 0.988, 0.661, 0.465, 0.883, 0.526, 0.213, 0.42, 0.993, 0.344, 0.547, 0.646, 0.117, 0.407, 0.011, 0.841, 0.521, 0.023, 0.17, 0.584, 0.367, 0.251, 0.434, 0.008, 0.434, 0.752, 0.7, 0.456, 0.106, 0.425, 0.532, 0.284, 0.83, 0.311, 0.234, 0.436, 0.917, 0.624, 0.054, 0.29, 0.639, 0.874, 0.204, 0.313, 0.035, 0.114, 0.064, 0.229, 0.52, 0.918, 0.917, 0.211, 0.995, 0.538, 0.294, 0.412, 0.21, 0.21, 0.629, 0.212, 0.01, 0.471, 0.564, 0.608, 0.701, 0.219, 0.927, 0.1, 0.6, 0.811, 0.468, 0.662, 0.791, 0.777, 0.837, 0.802, 0.024, 0.164, 0.002, 0.286, 0.544, 0.57, 0.439, 0.612, 0.102, 0.49, 0.921, 0.029, 0.942, 0.49, 0.266, 0.776, 0.822, 0.531, 0.583, 0.648, 0.218, 0.507, 0.079, 0.759, 0.277, 0.054, 0.762, 0.673, 0.494, 0.535, 0.211, 0.202, 0.424, 0.669, 0.485, 0.832, 0.004, 0.816, 0.733, 0.922, 0.519, 0.483, 0.442, 0.673, 0.09, 0.489, 0.148, 0.404, 0.194, 0.626, 0.024, 0.806, 0.256, 0.399, 0.381, 0.584, 0.021, 0.827, 0.237, 0.568, 0.609, 0.7, 0.908, 0.594, 0.458, 0.368, 0.392, 0.525, 0.382, 0.076, 0.865, 0.647, 0.614, 0.571, 0.077, 0.522, 0.019, 0.862, 0.111, 0.428, 0.649, 0.531, 0.559, 0.223, 0.2, 0.978, 0.245, 0.492, 0.174, 0.35, 0.291, 0.118, 0.543, 0.748, 0.22, 0.402, 0.805, 0.726, 0.247]
global q = [0.401, 0.982, 0.976, 0.04, 0.932, 0.935, 0.932, 0.905, 0.734, 0.949, 0.877, 0.964, 0.76, 0.643, 0.748, 0.932, 0.764, 0.812, 0.822, 0.097, 0.888, 0.994, 0.47, 0.999, 0.8, 0.641, 0.962, 0.644, 0.672, 0.494, 0.547, 0.638, 0.846, 0.887, 0.931, 0.281, 0.656, 0.494, 0.758, 0.594, 0.497, 0.996, 0.768, 0.91, 0.16, 0.977, 0.826, 0.861, 0.606, 0.998, 0.523, 0.995, 0.647, 0.99, 0.355, 0.712, 0.948, 0.962, 0.779, 0.765, 0.927, 0.779, 0.772, 0.836, 0.928, 0.891, 0.767, 0.569, 0.996, 0.913, 0.346, 0.668, 0.901, 0.831, 0.763, 0.547, 0.884, 0.972, 0.999, 0.922, 0.438, 0.573, 0.478, 0.985, 0.911, 0.962, 0.703, 0.926, 0.394, 0.785, 0.315, 0.042, 0.911, 0.967, 0.402, 0.994, 0.709, 0.932, 0.97, 0.395, 0.461, 0.787, 0.999, 0.859, 0.968, 0.751, 0.853, 0.709, 0.989, 0.824, 0.972, 0.611, 0.644, 0.84, 0.956, 0.904, 0.184, 0.663, 0.475, 0.259, 0.706, 0.517, 0.551, 0.855, 0.762, 0.876, 0.751, 0.537, 0.952, 0.564, 0.612, 0.938, 0.97, 0.216, 0.75, 0.689, 0.797, 0.554, 0.73, 0.673, 0.667, 0.924, 0.861, 0.975, 0.197, 0.578, 0.949, 0.966, 0.926, 0.837, 0.974, 0.805, 0.764, 0.405, 0.988, 0.822, 0.896, 0.958, 0.809, 0.838, 0.517, 0.993, 0.96, 0.63, 0.785, 0.384, 0.481, 0.484, 0.857, 0.678, 0.697, 0.739, 0.934, 0.523, 0.637, 0.682, 0.609, 0.775, 0.84, 0.763, 0.8, 0.308, 0.882, 0.799, 0.73, 0.926, 0.661, 0.86, 0.664, 0.993, 0.851, 0.64, 0.786, 0.719, 0.975, 0.443, 0.645, 0.22, 0.836, 0.227, 0.471, 0.529, 0.954, 0.964, 0.854, 0.999, 0.731, 0.418, 0.851, 0.514, 0.722, 0.671, 0.62, 0.655, 0.661, 0.919, 0.629, 0.807, 0.805, 0.989, 0.389, 0.62, 0.836, 0.9, 0.85, 0.896, 0.809, 0.862, 0.836, 0.639, 0.542, 0.899, 0.646, 0.941, 0.673, 0.529, 0.989, 0.564, 0.928, 0.982, 0.322, 0.943, 0.731, 0.99, 0.887, 0.85, 0.725, 0.647, 0.956, 0.949, 0.593, 0.861, 0.947, 0.999, 0.158, 0.855, 0.989, 0.597, 0.629, 0.836, 0.64, 0.902, 0.793, 0.78, 0.988, 0.342, 0.979, 0.983, 0.926, 0.666, 0.603, 0.532, 0.761, 0.616, 0.608, 0.989, 0.541, 0.346, 0.994, 0.693, 0.901, 0.475, 0.729, 0.432, 0.623, 0.754, 0.921, 0.96, 0.882, 0.696, 0.799, 0.952, 0.989, 0.526, 0.557, 0.551, 0.537, 0.637, 0.963, 0.95, 0.781, 0.677, 0.643, 0.688, 0.771, 0.914, 0.877, 0.827, 0.512, 0.718, 0.736, 0.867, 0.903, 0.833, 0.983, 0.399, 0.833, 0.534, 0.968, 0.912, 0.369, 0.744, 0.856, 0.761, 0.937, 0.917, 0.989, 0.802]
global origin = 1
global destination = 60