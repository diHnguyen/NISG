global arcs = [1 7; 1 11; 1 13; 1 23; 1 32; 1 33; 1 39; 1 52; 1 56; 2 17; 2 27; 2 36; 2 54; 2 55; 2 57; 3 2; 3 12; 3 14; 3 16; 3 20; 3 27; 3 31; 3 55; 3 60; 4 21; 4 25; 4 36; 5 4; 5 26; 5 33; 5 43; 5 48; 5 50; 6 8; 6 23; 6 46; 6 49; 6 55; 7 14; 7 48; 7 54; 7 56; 7 57; 8 5; 8 9; 8 24; 8 44; 8 46; 8 56; 8 60; 9 3; 9 7; 9 11; 9 17; 9 28; 9 36; 9 37; 9 38; 9 49; 10 2; 10 6; 10 12; 10 19; 10 43; 10 48; 11 10; 11 20; 11 25; 11 33; 11 36; 11 37; 11 40; 11 42; 11 50; 11 55; 12 15; 12 23; 12 32; 12 47; 12 49; 12 59; 13 4; 13 10; 13 17; 13 21; 13 24; 13 33; 13 40; 13 42; 14 20; 14 22; 14 26; 14 37; 14 47; 14 54; 14 57; 14 60; 15 4; 15 39; 15 44; 15 49; 16 5; 16 11; 16 18; 16 26; 16 34; 16 35; 16 54; 17 3; 17 22; 17 26; 17 32; 17 48; 17 54; 17 59; 18 7; 18 45; 18 49; 18 52; 19 23; 19 41; 19 57; 20 6; 20 7; 20 9; 20 45; 21 7; 21 19; 21 30; 21 32; 21 57; 22 27; 22 46; 23 5; 23 9; 23 28; 23 33; 24 16; 24 27; 24 38; 24 41; 24 42; 24 51; 24 54; 24 59; 25 23; 25 35; 25 46; 26 22; 27 5; 27 44; 28 20; 28 24; 28 40; 28 43; 28 48; 29 2; 29 3; 29 52; 30 17; 30 24; 30 45; 31 6; 31 9; 31 17; 31 28; 31 34; 31 44; 31 50; 32 11; 32 21; 32 23; 32 44; 32 52; 33 3; 33 10; 33 12; 33 45; 34 8; 34 13; 34 15; 34 24; 34 37; 34 40; 34 48; 34 49; 34 54; 35 34; 35 36; 35 38; 35 50; 36 10; 36 20; 36 26; 36 29; 36 46; 37 2; 37 7; 37 16; 37 32; 38 44; 38 52; 38 58; 39 10; 39 17; 39 25; 39 37; 39 40; 39 47; 39 55; 40 4; 40 20; 40 38; 40 45; 41 4; 41 30; 41 50; 41 51; 41 60; 42 12; 42 27; 42 28; 42 40; 43 10; 43 18; 43 32; 43 34; 43 37; 43 40; 43 47; 43 51; 43 58; 43 59; 44 5; 44 13; 44 23; 44 24; 44 34; 44 48; 44 58; 44 59; 45 7; 45 21; 45 33; 45 40; 45 46; 45 54; 45 60; 46 12; 46 13; 46 32; 46 47; 46 52; 47 3; 47 10; 47 13; 47 20; 47 53; 47 54; 47 58; 48 2; 48 11; 48 15; 48 31; 48 37; 48 45; 48 54; 48 58; 48 60; 49 4; 49 9; 49 11; 49 55; 49 59; 50 29; 50 52; 50 53; 51 32; 51 48; 51 50; 51 54; 52 12; 52 29; 52 32; 52 37; 52 45; 52 46; 52 47; 52 51; 52 55; 52 56; 53 12; 53 21; 54 16; 54 19; 54 26; 54 35; 54 36; 54 37; 54 52; 54 55; 54 60; 55 9; 55 13; 55 21; 55 23; 55 26; 55 48; 55 60; 56 7; 56 46; 56 53; 56 54; 57 12; 57 16; 57 24; 57 38; 57 52; 58 4; 58 8; 58 10; 58 26; 59 21; 59 22; 59 42; 59 47; 59 53; 59 58]
global d_x = [6.0, 10.0, 4.0, 10.0, 2.0, 5.0, 9.0, 8.0, 4.0, 3.0, 3.0, 4.0, 3.0, 7.0, 8.0, 5.0, 4.0, 2.0, 7.0, 5.0, 9.0, 6.0, 5.0, 7.0, 3.0, 10.0, 8.0, 6.0, 10.0, 6.0, 2.0, 9.0, 3.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 4.0, 4.0, 2.0, 2.0, 1.0, 7.0, 8.0, 7.0, 8.0, 7.0, 6.0, 5.0, 6.0, 10.0, 1.0, 4.0, 6.0, 5.0, 8.0, 9.0, 2.0, 3.0, 7.0, 1.0, 6.0, 3.0, 7.0, 3.0, 3.0, 7.0, 5.0, 3.0, 3.0, 4.0, 2.0, 6.0, 7.0, 6.0, 10.0, 6.0, 10.0, 10.0, 5.0, 5.0, 3.0, 7.0, 6.0, 7.0, 8.0, 5.0, 7.0, 6.0, 4.0, 7.0, 4.0, 1.0, 7.0, 10.0, 10.0, 9.0, 7.0, 3.0, 8.0, 9.0, 3.0, 6.0, 5.0, 10.0, 9.0, 9.0, 5.0, 6.0, 8.0, 4.0, 9.0, 2.0, 9.0, 6.0, 10.0, 1.0, 8.0, 3.0, 10.0, 10.0, 8.0, 7.0, 1.0, 5.0, 3.0, 4.0, 10.0, 8.0, 5.0, 3.0, 4.0, 7.0, 5.0, 3.0, 4.0, 3.0, 9.0, 6.0, 6.0, 2.0, 4.0, 5.0, 9.0, 3.0, 1.0, 4.0, 7.0, 9.0, 8.0, 3.0, 8.0, 7.0, 3.0, 1.0, 3.0, 8.0, 10.0, 3.0, 7.0, 1.0, 9.0, 8.0, 3.0, 8.0, 8.0, 4.0, 7.0, 2.0, 10.0, 2.0, 10.0, 10.0, 7.0, 8.0, 2.0, 1.0, 2.0, 4.0, 5.0, 10.0, 8.0, 9.0, 9.0, 4.0, 4.0, 7.0, 1.0, 6.0, 10.0, 10.0, 8.0, 4.0, 3.0, 8.0, 4.0, 9.0, 2.0, 4.0, 3.0, 6.0, 5.0, 6.0, 3.0, 3.0, 9.0, 1.0, 9.0, 4.0, 10.0, 10.0, 10.0, 3.0, 3.0, 6.0, 9.0, 1.0, 3.0, 9.0, 5.0, 3.0, 8.0, 8.0, 7.0, 5.0, 2.0, 10.0, 1.0, 4.0, 8.0, 10.0, 2.0, 2.0, 5.0, 5.0, 5.0, 10.0, 2.0, 4.0, 7.0, 3.0, 9.0, 4.0, 1.0, 8.0, 8.0, 6.0, 10.0, 7.0, 4.0, 3.0, 7.0, 8.0, 2.0, 1.0, 4.0, 5.0, 6.0, 10.0, 8.0, 9.0, 4.0, 2.0, 4.0, 3.0, 7.0, 3.0, 10.0, 4.0, 5.0, 2.0, 2.0, 3.0, 5.0, 4.0, 10.0, 4.0, 3.0, 7.0, 4.0, 2.0, 9.0, 10.0, 3.0, 4.0, 3.0, 5.0, 8.0, 4.0, 10.0, 4.0, 7.0, 4.0, 6.0, 2.0, 5.0, 4.0, 1.0, 9.0, 1.0, 4.0, 9.0, 2.0, 6.0, 5.0, 9.0, 1.0, 7.0, 8.0, 8.0, 4.0, 5.0, 6.0, 1.0, 3.0, 7.0, 4.0, 4.0, 5.0, 5.0, 8.0, 9.0, 1.0, 1.0, 4.0, 6.0]
global b_x = 5
global d_y = [6.0, 1.0, 7.0, 5.0, 3.0, 7.0, 4.0, 10.0, 6.0, 2.0, 4.0, 1.0, 4.0, 2.0, 10.0, 1.0, 5.0, 7.0, 7.0, 3.0, 4.0, 4.0, 7.0, 5.0, 2.0, 10.0, 8.0, 6.0, 1.0, 10.0, 9.0, 4.0, 3.0, 10.0, 5.0, 3.0, 7.0, 6.0, 3.0, 9.0, 3.0, 7.0, 3.0, 6.0, 4.0, 9.0, 3.0, 8.0, 6.0, 8.0, 9.0, 3.0, 4.0, 4.0, 2.0, 2.0, 7.0, 1.0, 8.0, 6.0, 10.0, 3.0, 6.0, 9.0, 1.0, 7.0, 9.0, 9.0, 5.0, 2.0, 7.0, 10.0, 7.0, 3.0, 6.0, 5.0, 6.0, 2.0, 7.0, 9.0, 3.0, 4.0, 1.0, 1.0, 3.0, 7.0, 10.0, 10.0, 5.0, 2.0, 7.0, 2.0, 5.0, 5.0, 1.0, 2.0, 7.0, 2.0, 3.0, 2.0, 1.0, 5.0, 10.0, 3.0, 6.0, 8.0, 7.0, 5.0, 2.0, 6.0, 9.0, 7.0, 9.0, 4.0, 3.0, 1.0, 9.0, 9.0, 2.0, 5.0, 7.0, 3.0, 10.0, 4.0, 4.0, 5.0, 4.0, 6.0, 5.0, 8.0, 4.0, 5.0, 8.0, 9.0, 8.0, 2.0, 1.0, 8.0, 8.0, 8.0, 5.0, 4.0, 7.0, 10.0, 6.0, 8.0, 10.0, 10.0, 3.0, 10.0, 6.0, 1.0, 6.0, 2.0, 3.0, 8.0, 2.0, 3.0, 1.0, 1.0, 7.0, 2.0, 6.0, 7.0, 5.0, 9.0, 3.0, 2.0, 3.0, 5.0, 9.0, 4.0, 6.0, 10.0, 7.0, 8.0, 2.0, 8.0, 1.0, 10.0, 2.0, 1.0, 9.0, 8.0, 1.0, 1.0, 6.0, 8.0, 8.0, 4.0, 1.0, 2.0, 4.0, 9.0, 2.0, 3.0, 9.0, 10.0, 1.0, 5.0, 1.0, 5.0, 8.0, 9.0, 8.0, 10.0, 3.0, 10.0, 8.0, 4.0, 10.0, 8.0, 10.0, 8.0, 10.0, 4.0, 10.0, 9.0, 9.0, 5.0, 10.0, 4.0, 8.0, 1.0, 9.0, 3.0, 4.0, 8.0, 1.0, 6.0, 4.0, 1.0, 3.0, 9.0, 5.0, 7.0, 7.0, 5.0, 7.0, 7.0, 4.0, 9.0, 6.0, 7.0, 8.0, 7.0, 10.0, 5.0, 5.0, 6.0, 9.0, 8.0, 4.0, 8.0, 4.0, 6.0, 7.0, 10.0, 4.0, 2.0, 6.0, 10.0, 10.0, 9.0, 7.0, 9.0, 10.0, 1.0, 9.0, 4.0, 3.0, 10.0, 9.0, 1.0, 10.0, 4.0, 5.0, 6.0, 8.0, 1.0, 7.0, 7.0, 8.0, 5.0, 5.0, 10.0, 10.0, 2.0, 7.0, 10.0, 1.0, 7.0, 9.0, 1.0, 5.0, 6.0, 5.0, 3.0, 10.0, 9.0, 4.0, 5.0, 4.0, 2.0, 1.0, 5.0, 9.0, 6.0, 6.0, 9.0, 10.0, 2.0, 6.0, 10.0, 2.0, 2.0, 9.0, 4.0, 3.0, 2.0, 2.0, 6.0, 10.0, 1.0, 10.0, 4.0, 1.0, 3.0]
global b_y = 10
global p = [0.421, 0.588, 0.402, 0.512, 0.378, 0.031, 0.094, 0.85, 0.793, 0.717, 0.499, 0.643, 0.471, 0.74, 0.219, 0.685, 0.714, 0.995, 0.188, 0.606, 0.622, 0.468, 0.744, 0.449, 0.232, 0.614, 0.45, 0.463, 0.866, 0.505, 0.219, 0.816, 0.757, 0.365, 0.128, 0.776, 0.788, 0.084, 0.541, 0.997, 0.453, 0.568, 0.84, 0.76, 0.57, 0.502, 0.594, 0.557, 0.091, 0.783, 0.99, 0.894, 0.986, 0.8, 0.561, 0.748, 0.752, 0.985, 0.97, 0.517, 0.883, 0.696, 0.37, 0.185, 0.096, 0.861, 0.417, 0.344, 0.209, 0.544, 0.311, 0.647, 0.185, 0.989, 0.305, 0.199, 0.322, 0.801, 0.379, 0.158, 0.58, 0.744, 0.44, 0.464, 0.285, 0.533, 0.548, 0.649, 0.157, 0.397, 0.528, 0.529, 0.552, 0.652, 0.17, 0.773, 0.291, 0.893, 0.306, 0.465, 0.551, 0.695, 0.009, 0.998, 0.142, 0.463, 0.775, 0.537, 0.624, 0.56, 0.912, 0.225, 0.969, 0.903, 0.193, 0.746, 0.133, 0.221, 0.351, 0.636, 0.95, 0.637, 0.314, 0.81, 0.657, 0.035, 0.619, 0.121, 0.667, 0.774, 0.123, 0.298, 0.923, 0.032, 0.638, 0.309, 0.645, 0.562, 0.877, 0.039, 0.722, 0.097, 0.151, 0.066, 0.543, 0.858, 0.363, 0.268, 0.477, 0.575, 0.448, 0.229, 0.358, 0.784, 0.12, 0.852, 0.268, 0.395, 0.102, 0.805, 0.319, 0.789, 0.972, 0.569, 0.429, 0.121, 0.424, 0.906, 0.682, 0.293, 0.577, 0.829, 0.039, 0.237, 0.171, 0.54, 0.959, 0.599, 0.494, 0.726, 0.702, 0.084, 0.243, 0.888, 0.877, 0.82, 0.22, 0.892, 0.057, 0.469, 0.466, 0.007, 0.732, 0.73, 0.51, 0.733, 0.058, 0.143, 0.046, 0.331, 0.049, 0.323, 0.055, 0.917, 0.363, 0.97, 0.916, 0.751, 0.521, 0.597, 0.699, 0.346, 0.352, 0.092, 0.967, 0.376, 0.669, 0.073, 0.335, 0.507, 0.672, 0.892, 0.829, 0.2, 0.927, 0.058, 0.17, 0.233, 0.979, 0.42, 0.735, 0.349, 0.66, 0.702, 0.263, 0.443, 0.733, 0.566, 0.581, 0.695, 0.132, 0.271, 0.333, 0.899, 0.914, 0.125, 0.536, 0.352, 0.48, 0.053, 0.51, 0.349, 0.381, 0.338, 0.309, 0.46, 0.94, 0.45, 0.057, 0.286, 0.208, 0.525, 0.003, 0.197, 0.71, 0.772, 0.544, 0.299, 0.958, 0.032, 0.544, 0.244, 0.381, 0.298, 0.781, 0.949, 0.036, 0.791, 0.905, 0.389, 0.944, 0.867, 0.04, 0.318, 0.964, 0.115, 0.475, 0.304, 0.643, 0.374, 0.559, 0.859, 0.139, 0.431, 0.75, 0.679, 0.698, 0.829, 0.902, 0.019, 0.372, 0.149, 0.41, 0.941, 0.706, 0.086, 0.916, 0.451, 0.741, 0.633, 0.967, 0.41, 0.658, 0.791, 0.86, 0.589, 0.279, 0.212, 0.284, 0.633, 0.726, 0.928, 0.662, 0.175, 0.768, 0.262, 0.634, 0.964]
global q = [0.921, 0.663, 0.695, 0.781, 0.939, 0.059, 0.978, 0.984, 0.959, 0.838, 0.789, 0.755, 0.546, 0.802, 0.242, 0.882, 0.74, 0.998, 0.925, 0.706, 0.779, 0.796, 0.849, 0.995, 0.564, 0.895, 0.788, 0.934, 0.882, 0.745, 0.993, 0.879, 0.835, 0.758, 0.517, 0.871, 0.899, 0.853, 0.758, 0.997, 0.852, 0.907, 0.9, 0.933, 0.879, 0.973, 0.672, 0.992, 0.201, 0.948, 0.992, 0.951, 0.995, 0.978, 0.6, 0.777, 0.752, 0.992, 0.971, 0.751, 0.93, 0.844, 0.74, 0.667, 0.669, 0.942, 0.448, 0.792, 0.639, 0.886, 0.529, 0.706, 0.66, 0.994, 0.923, 0.708, 0.439, 0.821, 0.841, 0.537, 0.872, 0.991, 0.521, 0.756, 0.371, 0.854, 0.608, 0.845, 0.295, 0.678, 0.908, 0.806, 0.784, 0.92, 0.508, 0.975, 0.34, 0.944, 0.367, 0.855, 0.603, 0.725, 0.394, 0.999, 0.582, 0.941, 0.958, 0.542, 0.992, 0.995, 0.992, 0.501, 0.972, 0.905, 0.652, 0.905, 0.544, 0.589, 0.637, 0.71, 0.992, 0.875, 0.941, 0.917, 0.857, 0.451, 0.894, 0.676, 0.696, 0.816, 0.302, 0.808, 0.982, 0.516, 0.9, 0.877, 0.671, 0.636, 0.883, 0.385, 0.997, 0.115, 0.439, 0.447, 0.547, 0.874, 0.583, 0.967, 0.557, 0.841, 0.525, 0.326, 0.586, 0.962, 0.17, 0.976, 0.777, 0.908, 0.667, 0.911, 0.796, 0.817, 0.98, 0.661, 0.59, 0.159, 0.526, 0.992, 0.802, 0.499, 0.902, 0.864, 0.654, 0.481, 0.932, 0.594, 0.961, 0.682, 0.723, 0.983, 0.764, 0.459, 0.492, 0.983, 0.909, 0.971, 0.742, 0.91, 0.446, 0.901, 0.7, 0.369, 0.876, 0.895, 0.772, 0.838, 0.494, 0.958, 0.446, 0.578, 0.422, 0.777, 0.224, 0.949, 0.558, 0.978, 0.964, 0.963, 0.68, 0.953, 0.968, 0.369, 0.88, 0.791, 0.996, 0.602, 0.784, 0.847, 0.93, 0.528, 0.88, 0.954, 0.951, 0.722, 0.972, 0.53, 0.449, 0.497, 0.98, 0.756, 0.755, 0.548, 0.874, 0.924, 0.557, 0.902, 0.964, 0.619, 0.611, 0.96, 0.827, 0.768, 0.48, 0.941, 0.95, 0.841, 0.646, 0.489, 0.934, 0.28, 0.585, 0.697, 0.601, 0.543, 0.497, 0.872, 0.949, 0.85, 0.689, 0.9, 0.386, 0.701, 0.306, 0.564, 0.786, 0.827, 0.866, 0.9, 0.964, 0.545, 0.569, 0.597, 0.914, 0.918, 0.962, 0.959, 0.5, 0.848, 0.984, 0.425, 0.962, 0.874, 0.666, 0.92, 0.986, 0.948, 0.98, 0.684, 0.646, 0.562, 0.974, 0.927, 0.235, 0.861, 0.853, 0.89, 0.999, 0.866, 0.993, 0.121, 0.48, 0.378, 0.767, 0.961, 0.815, 0.887, 0.949, 0.919, 0.774, 0.806, 0.999, 0.683, 0.999, 0.999, 0.986, 0.653, 0.763, 0.721, 0.449, 0.738, 0.882, 0.967, 0.924, 0.902, 0.993, 0.76, 0.826, 0.997]
global origin = 1
global destination = 60